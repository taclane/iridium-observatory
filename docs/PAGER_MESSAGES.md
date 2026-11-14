# Pager (MSG) traffic: what the messages are

The simplex messaging channel (`MSG` frames, 1626.0–1626.5 MHz) carries one-way
pages. `frame_decode.c` decodes them as format 5 (7-bit ASCII) or format 3
(BCD), and the archive keeps them as `{"e":"pager", ric, fmt, ok, txt}` events.
This document records what that traffic turned out to be. It is based on one
station's archive from 2026-11-12 to 2026-11-14: 5,827 pager events, 63 RICs,
all format 5.

## Two kinds of traffic

About a quarter of the events are **plain text**, sent to a few low RICs
(93826, 100524 and a handful of others):

| Pattern | Notes |
|---|---|
| `nm@smtp.iwcs.net Test message sent from SNOC-CRON -> 8816310xx` | Scheduled test page from the network operations centre (SNOC) |
| `NNN@msg.iridium.com on <date> UTC <year>` | Hourly timestamped heartbeat at hh:15 |
| `EMSS Indoor Simplex Test` | Test page |
| Occasional operator emails | Person-to-pager email; do not republish |

Short versions of the same texts (`nm@sm`, `646@msg.iridium.com on`) are
receptions that were cut off partway through.

The rest is **base64 fragments of larger binary messages**. They look like
`CgAMCSZviS9zf*JeNu*J…` and go to a pool of RICs in the 3375xxx/352xxxx/98xxx
ranges. Text of this kind made up about 70% of the pager events.

## Fragment format

The pager alphabet uses two substitutions: `*` stands for base64 `/` and `-`
for the `=` pad. With those mapped back, each fragment decodes to:

```
byte 0-1   message id, little-endian (byte 1 is almost always 0)
byte 2     total fragments in the message
byte 3     fragment index, 0 .. total-1
byte 4     payload length: 38 for a full fragment (60 base64 chars), less for the last
5..5+len   payload
last 2     CRC-16/XMODEM (poly 0x1021, init 0), big-endian, over bytes 0 .. 5+len-1
```

The CRC confirms this layout: 153 of the 158 fragments whose length byte
matched their size passed it. The rest are consistent with bit errors.

## Transmission pattern

- **Each fragment goes out twice**, about 3 minutes apart (median gap 181–196 s).
- **The fragments of one message are spread over several RICs.** The schedule
  works through them in turn, about 9 minutes per RIC. One fragment can reach
  3–4 RICs in a cycle, and the fragment-to-RIC mapping changes from cycle to
  cycle. The RICs look like delivery slots that a group of terminals listens on,
  not individual recipients.
- **Whole messages are rebroadcast unchanged on a fixed cycle**, mostly about
  4 h (message 5 came back at exactly 4.00 h intervals) and sometimes about 21 h.
  The fragments are byte-identical each time.

Taken together, this is a looping broadcast to a fleet of terminals: a carousel
of configuration, keys, schedules or a bulletin. It is not a stream of new
data. "The same text again tens of minutes later" is the carousel moving on to
the next RIC, or starting its next cycle.

## Payload

Reassembled messages so far:

| id | fragments | bytes | first two bytes |
|---|---|---|---|
| 5 | 9 | 308 | `02 82` |
| 10 | 12 | 429 | `02 83` |
| 11 | 8 | 288 | `02 82` |
| 12 | 9 | 306 | `02 82` |
| 314 | 8 | 292 | `02 82` |

There are also eight messages of 53–54 fragments (about 2 KB each, ids 3, 4,
61, 96, 102, 110, 116, 207). None of them were complete yet.

- Every complete payload starts with `02 82` or `02 83`, followed by a byte that
  varies. This is probably a small plaintext header (a version or type byte and
  flags) in front of the body.
- After that the data cannot be told apart from random. Entropy is 7.22–7.53
  bits/byte; uniformly random 300-byte strings average 7.29.
- Lengths are not multiples of 8 or 16, which rules out a padded block cipher.
  A stream or AEAD mode (CTR, GCM), or compressed data, fits.
- Repeating identical ciphertext is normal for a broadcast that was encrypted
  once. It says nothing about how weak the encryption is.

Conclusion so far: the payload is encrypted or compressed application data, and
nothing more can be read from it. Its framing and broadcast schedule are fully
understood.

## Tooling

`observatory/tools/pager_reasm.py` reads the archive, keeps the fragments that
pass the CRC, and groups them by (id, total). It reports completeness, broadcast
cycles and the RICs used. With `--dump DIR` it writes each complete payload to a
`.bin` file.

```
python3 observatory/tools/pager_reasm.py --days 3 --dump /tmp/pager-msgs
```

## Decoder note

Until this branch, the archived `ok` flag was 0 for every pager message.
`parse_msg` ran the packet checksum over all message blocks, but iridium-toolkit
(`msg_checksum(self.blocks[1:])`) excludes the first block, which holds the
start of the RIC. Archives written before the fix have `ok:0` throughout;
treat the flag as meaningful only from the fix onward.
