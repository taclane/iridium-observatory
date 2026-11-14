#!/usr/bin/env python3
"""
pager_reasm.py - reassemble the fragmented binary pager broadcasts in the archive.

Much of the simplex pager traffic is not human text: it is base64 fragments of
larger binary messages. Each fragment, once '*' -> '/' and '-' -> '=' are
mapped back to standard base64, decodes to:

    byte 0-1   message id (little-endian; byte 1 is almost always 0)
    byte 2     total fragments in the message
    byte 3     fragment index (0 .. total-1)
    byte 4     payload length (38 for a full fragment, less for the last)
    5..5+len   payload
    last 2     CRC-16/XMODEM, big-endian, over bytes 0 .. 5+len-1

The fragments of one message are spread over a pool of RICs, so the RIC is not
part of the key. Each fragment is sent twice about 3 minutes apart, and a whole
message is rebroadcast unchanged on a multi-hour cycle. This tool groups the
CRC-valid fragments by (id, total), reports how complete each message is and
how often it was broadcast, and can write the complete payloads out.

    python3 observatory/tools/pager_reasm.py
    python3 observatory/tools/pager_reasm.py --archive /var/lib/iridium-observatory/archive \\
        --days 3 --dump /tmp/pager-msgs

Copyright (c) 2026 Mike Brown
SPDX-License-Identifier: GPL-3.0-or-later
"""

import argparse
import base64
import binascii
import collections
import glob
import hashlib
import json
import math
import os
import time

EPISODE_GAP_S = 30 * 60   # broadcasts of one fragment further apart than this are separate cycles


def decode_fragment(txt):
    """Return (id, total, index, payload) for a CRC-valid fragment, else None."""
    if not txt or " " in txt or "@" in txt:
        return None
    s = txt.replace("*", "/").replace("-", "=")
    try:
        b = base64.b64decode(s + "=" * (-len(s) % 4), validate=True)
    except (binascii.Error, ValueError):
        return None
    if len(b) < 7 or b[4] + 7 != len(b):
        return None
    if binascii.crc_hqx(b[:-2], 0) != int.from_bytes(b[-2:], "big"):
        return None
    mid, total, idx = b[0] | (b[1] << 8), b[2], b[3]
    if total == 0 or idx >= total:
        return None
    return mid, total, idx, bytes(b[5:-2])


def entropy(data):
    if not data:
        return 0.0
    c = collections.Counter(data)
    n = len(data)
    return -sum(v / n * math.log2(v / n) for v in c.values())


def episodes(times):
    times = sorted(times)
    starts = [times[0]]
    for a, b in zip(times, times[1:]):
        if b - a > EPISODE_GAP_S:
            starts.append(b)
    return starts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--archive", default="/var/lib/iridium-observatory/archive")
    ap.add_argument("--days", type=int, default=0, help="only the newest N day files (0 = all)")
    ap.add_argument("--dump", metavar="DIR", help="write each complete message payload to DIR")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.archive, "iridium-*.jsonl")))
    if args.days:
        files = files[-args.days:]

    stats = collections.Counter()
    # (id, total) -> index -> payload -> list of (t, ric)
    msgs = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    for fn in files:
        with open(fn, errors="replace") as f:
            for line in f:
                if '"pager"' not in line:
                    continue
                try:
                    ev = json.loads(line)
                except ValueError:
                    continue
                if ev.get("e") != "pager":
                    continue
                stats["pager events"] += 1
                txt = ev.get("txt", "")
                if " " in txt or "@" in txt:
                    stats["plain text"] += 1
                    continue
                frag = decode_fragment(txt)
                if frag is None:
                    stats["rejected (not base64 / short / CRC)"] += 1
                    continue
                stats["valid fragments"] += 1
                mid, total, idx, payload = frag
                msgs[(mid, total)][idx][payload].append((ev["t"] / 1000.0, ev.get("ric")))

    for k, v in stats.items():
        print(f"{k:>40}: {v}")
    print()

    if args.dump:
        os.makedirs(args.dump, exist_ok=True)

    hdr = f"{'id':>5} {'frags':>7} {'bytes':>6} {'cycles':>6} {'rics':>4}  {'first seen (UTC)':16}  {'entropy':>7}  sha1/notes"
    print(hdr)
    print("-" * len(hdr))
    complete = 0
    for (mid, total), frags in sorted(msgs.items(), key=lambda kv: min(t for p in kv[1].values() for h in p.values() for t, _ in h)):
        all_hits = [h for p in frags.values() for hits in p.values() for h in hits]
        first = min(t for t, _ in all_hits)
        rics = {r for _, r in all_hits}
        # Cycle count from the fragment seen most often, so one missed fragment
        # does not undercount.
        cycles = max(len(episodes([t for t, _ in hits])) for p in frags.values() for hits in p.values())
        conflicts = sorted(i for i, p in frags.items() if len(p) > 1)
        have = len(frags)
        notes = []
        if conflicts:
            notes.append(f"id reused? differing payloads at idx {conflicts[:6]}")
        line = f"{mid:>5} {have:>3}/{total:<3} "
        if have == total and not conflicts:
            complete += 1
            data = b"".join(next(iter(frags[i])) for i in range(total))
            line += f"{len(data):>6} {cycles:>6} {len(rics):>4}  {time.strftime('%m-%d %H:%M', time.gmtime(first)):16}  {entropy(data):7.3f}  {hashlib.sha1(data).hexdigest()[:12]}"
            if args.dump:
                path = os.path.join(args.dump, f"msg-{mid:04x}-{total}-{hashlib.sha1(data).hexdigest()[:12]}.bin")
                with open(path, "wb") as f:
                    f.write(data)
        else:
            missing = [i for i in range(total) if i not in frags]
            if missing:
                notes.insert(0, f"missing {missing[:10]}{'...' if len(missing) > 10 else ''}")
            line += f"{'-':>6} {cycles:>6} {len(rics):>4}  {time.strftime('%m-%d %H:%M', time.gmtime(first)):16}  {'-':>7}"
        print(line + ("  " + "; ".join(notes) if notes else ""))

    print(f"\n{complete} of {len(msgs)} messages complete")


if __name__ == "__main__":
    main()
