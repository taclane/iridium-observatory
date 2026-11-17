/*
 * Overlap-save burst channelizer (ka9q-radio style fast convolution)
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * See os_chan.h for what this is for. Notation:
 *   D  decimation, in_rate / out_rate (40 at 10 MS/s -> 250 kS/s)
 *   M  output bins per block, OS_CHAN_M (1024: 244 Hz bins)
 *   N  block FFT size, M * D (40960 at 10 MS/s)
 *   L  new input samples per block, 15N/16; V = N - L overlap
 *
 * Block b transforms input [bL - V, bL + L). With the filter's taps h[0..P)
 * zero-padded to N, the circular convolution of that block with h is the
 * linear one for the last L positions (V >= P - 1), and decimating it by D
 * is an M-point inverse FFT of the M bins around the burst's frequency:
 *     z[j] = sum_t h[t] x[B + jD + tau - t] e^{-2 pi i c (jD + tau - t) / N}
 * for block start B = bL - V, bin shift c and a time advance tau < D (a
 * phase ramp on the bins). The old path's out[k] is that same sum at
 * m = start + kD + P - 1, phase-referenced to start instead of B, so tau is
 * chosen to put the D-sample grid on m and each block's output is turned by
 * e^{2 pi i c (start - B) / N}. h is symmetric, so correlation (the old FIR's
 * form) and convolution are the same filter.
 */

#define _GNU_SOURCE             /* pthread_setname_np */
#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "os_chan.h"
#include "fir_filter.h"
#include "fftw_lock.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define OS_CHAN_M      1024
/* Cached block FFTs: 64 blocks of 3.84 ms (at 10 MS/s) span ~245 ms, well
 * past how far apart in time the bursts of one emission batch start. A miss
 * costs only a recompute. */
#define OS_CHAN_SLOTS  64

/* With the block thread on, how many blocks past the detector's position it
 * may transform. A burst reaches back at most the ringbuffer (about 33
 * blocks at 10 MS/s), so AHEAD plus that stays under OS_CHAN_SLOTS and the
 * thread never reuses a slot a burst is still reading; get_block checks. */
#define OS_CHAN_AHEAD  8
/* ...and how far behind it the thread stops bothering: a block that old is
 * one no burst still to come will ask for (when the detector outruns the
 * thread, reading a file faster than real time). */
#define OS_CHAN_BEHIND (OS_CHAN_SLOTS - OS_CHAN_AHEAD - 8)

struct os_chan {
    int D, M, N, L, V, ntaps;
    uint64_t delay_ns;

    /* G[tau * M + idx]: filter response at bin idx (k' = idx - M if
     * idx >= M/2), times the tau-sample advance ramp, times 1/N. */
    float complex *G;

    fftwf_plan fwd;             /* N-point, fin -> a slot or fout */
    float complex *fin, *fout;
    fftwf_plan inv;             /* M-point, zin -> zout */
    float complex *zin, *zout;

    float complex *slots;       /* OS_CHAN_SLOTS * N */
    int64_t tag[OS_CHAN_SLOTS]; /* block index held, -1 if none */

    os_chan_read_fn read;
    void *read_ctx;

    uint64_t n_fft, n_hit;

    /* Block thread (os_chan_start_thread). Block FFTs are most of the
     * detector thread's work, and nearly every block is needed, so a second
     * thread transforms each one as soon as its input is in; get_block then
     * finds it cached. Everything below is guarded by mu. */
    int threaded, stop;
    pthread_t th;
    pthread_mutex_t mu;
    pthread_cond_t cv;          /* any change below: one condition, broadcast */
    unsigned char busy[OS_CHAN_SLOTS];  /* slot being transformed (tag says which) */
    int64_t pos_b;              /* block under the detector's state machine */
    uint64_t have_hi;           /* input received: [0, have_hi) */
    uint64_t keep_lo;           /* the detector may overwrite input below this */
    int64_t next_b;             /* next block for the thread */
    int reading;                /* thread copying input [read_lo, ...) */
    uint64_t read_lo;
    float complex *tfin;        /* the thread's own input buffer */
};

os_chan_t *os_chan_create(int in_rate, int out_rate,
                          os_chan_read_fn read, void *read_ctx) {
    if (in_rate <= 0 || out_rate <= 0 || in_rate % out_rate != 0)
        return NULL;
    int D = in_rate / out_rate;
    if (D < 2)
        return NULL;

    os_chan_t *c = calloc(1, sizeof(*c));
    c->D = D;
    c->M = OS_CHAN_M;
    c->N = c->M * D;
    /* The overlap only has to hold the filter (V >= ntaps - 1: 2560 against
     * 801 at 10 MS/s, and 64D against ~20D at any rate), and every sample of
     * overlap is transformed twice, so keep it small. */
    c->V = c->N / 16;
    c->L = c->N - c->V;
    c->read = read;
    c->read_ctx = read_ctx;
    for (int i = 0; i < OS_CHAN_SLOTS; i++)
        c->tag[i] = -1;

    /* The same filter the time-domain path uses for a one-stage decimation. */
    int ntaps;
    float *taps = lpf_taps(&ntaps, 1.0f, (float)in_rate,
                           out_rate * 0.4f, out_rate * 0.2f);
    c->ntaps = ntaps;
    c->delay_ns = (uint64_t)((double)(ntaps / 2) / in_rate * 1e9);
    if (ntaps - 1 > c->V) {         /* can't happen for these designs */
        free(taps);
        free(c);
        return NULL;
    }

    c->fin = fftwf_alloc_complex(c->N);
    c->fout = fftwf_alloc_complex(c->N);
    c->slots = fftwf_alloc_complex((size_t)OS_CHAN_SLOTS * c->N);
    c->zin = fftwf_alloc_complex(c->M);
    c->zout = fftwf_alloc_complex(c->M);
    c->G = malloc(sizeof(float complex) * (size_t)D * c->M);

    fftw_lock();
    c->fwd = fftwf_plan_dft_1d(c->N, (fftwf_complex *)c->fin,
                               (fftwf_complex *)c->fout,
                               FFTW_FORWARD, FFTW_MEASURE);
    c->inv = fftwf_plan_dft_1d(c->M, (fftwf_complex *)c->zin,
                               (fftwf_complex *)c->zout,
                               FFTW_BACKWARD, FFTW_MEASURE);
    fftw_unlock();

    /* Filter response at the M bins around 0 (planning may have scribbled
     * on fin, so fill it after). */
    memset(c->fin, 0, sizeof(float complex) * c->N);
    for (int t = 0; t < ntaps; t++)
        c->fin[t] = taps[t];
    free(taps);
    fftwf_execute(c->fwd);
    for (int tau = 0; tau < D; tau++) {
        for (int idx = 0; idx < c->M; idx++) {
            int kp = idx < c->M / 2 ? idx : idx - c->M;
            double complex ramp = cexp(I * 2.0 * M_PI * kp * tau / c->N);
            c->G[(size_t)tau * c->M + idx] =
                (float complex)(c->fout[(kp + c->N) % c->N] * ramp / c->N);
        }
    }
    return c;
}

void os_chan_destroy(os_chan_t *c) {
    if (!c) return;
    if (c->threaded) {
        pthread_mutex_lock(&c->mu);
        c->stop = 1;
        pthread_cond_broadcast(&c->cv);
        pthread_mutex_unlock(&c->mu);
        pthread_join(c->th, NULL);
        pthread_mutex_destroy(&c->mu);
        pthread_cond_destroy(&c->cv);
        fftwf_free(c->tfin);
    }
    fftw_lock();
    fftwf_destroy_plan(c->fwd);
    fftwf_destroy_plan(c->inv);
    fftw_unlock();
    fftwf_free(c->fin);
    fftwf_free(c->fout);
    fftwf_free(c->slots);
    fftwf_free(c->zin);
    fftwf_free(c->zout);
    free(c->G);
    free(c);
}

int os_chan_bins_per(const os_chan_t *c, int fft_size) {
    return fft_size > 0 && c->N % fft_size == 0 ? c->N / fft_size : 0;
}

int os_chan_fft_size(const os_chan_t *c) {
    return c->N;
}

void os_chan_stats(const os_chan_t *c, uint64_t *ffts, uint64_t *hits) {
    os_chan_t *m = (os_chan_t *)c;      /* the lock, not the counters */
    if (c->threaded)
        pthread_mutex_lock(&m->mu);
    *ffts = c->n_fft;
    *hits = c->n_hit;
    if (c->threaded)
        pthread_mutex_unlock(&m->mu);
}

/* FFT of block b, from the cache if it's there. A block whose input wasn't
 * all available (start of stream, aged out, not received yet) is computed
 * into fout and not cached, so a later burst can't get the zero-filled one. */
static const float complex *get_block_threaded(os_chan_t *c, int64_t b);

static const float complex *get_block(os_chan_t *c, int64_t b) {
    if (c->threaded)
        return get_block_threaded(c, b);
    int slot = (int)(b % OS_CHAN_SLOTS);
    float complex *dst = &c->slots[(size_t)slot * c->N];
    if (c->tag[slot] == b) {
        c->n_hit++;
        return dst;
    }

    int64_t B = b * c->L - c->V;
    int complete;
    if (B >= 0) {
        complete = c->read(c->read_ctx, (uint64_t)B, (size_t)c->N, c->fin);
    } else {
        size_t skip = (size_t)(-B);
        memset(c->fin, 0, sizeof(float complex) * skip);
        c->read(c->read_ctx, 0, (size_t)c->N - skip, c->fin + skip);
        complete = 0;
    }

    if (!complete)
        dst = c->fout;
    fftwf_execute_dft(c->fwd, (fftwf_complex *)c->fin, (fftwf_complex *)dst);
    c->n_fft++;
    if (complete)
        c->tag[slot] = b;
    return dst;
}

/* The detector thread's get_block when the block thread runs. A slot the
 * thread is still transforming is waited for; a block it hasn't reached is
 * claimed and transformed here, as without the thread. A block so far
 * behind that the thread may reuse its slot while the caller reads it is
 * transformed into fout and not cached (it shouldn't happen: see AHEAD). */
static const float complex *get_block_threaded(os_chan_t *c, int64_t b) {
    int slot = (int)(b % OS_CHAN_SLOTS);
    float complex *dst = &c->slots[(size_t)slot * c->N];
    pthread_mutex_lock(&c->mu);
    int cacheable = b + OS_CHAN_SLOTS > c->pos_b + OS_CHAN_AHEAD;
    for (;;) {
        if (!cacheable)
            break;
        if (c->tag[slot] == b && !c->busy[slot]) {
            c->n_hit++;
            pthread_mutex_unlock(&c->mu);
            return dst;
        }
        if (!c->busy[slot])
            break;
        pthread_cond_wait(&c->cv, &c->mu);     /* the thread has this slot */
    }
    if (cacheable) {
        c->tag[slot] = b;
        c->busy[slot] = 1;
    }
    pthread_mutex_unlock(&c->mu);

    int64_t B = b * c->L - c->V;
    int complete;
    if (B >= 0) {
        complete = c->read(c->read_ctx, (uint64_t)B, (size_t)c->N, c->fin);
    } else {
        size_t skip = (size_t)(-B);
        memset(c->fin, 0, sizeof(float complex) * skip);
        c->read(c->read_ctx, 0, (size_t)c->N - skip, c->fin + skip);
        complete = 0;
    }
    if (!complete || !cacheable)
        dst = c->fout;
    fftwf_execute_dft(c->fwd, (fftwf_complex *)c->fin, (fftwf_complex *)dst);

    pthread_mutex_lock(&c->mu);
    c->n_fft++;
    if (cacheable) {
        c->busy[slot] = 0;
        if (!complete)
            c->tag[slot] = -1;
        pthread_cond_broadcast(&c->cv);
    }
    pthread_mutex_unlock(&c->mu);
    return dst;
}

static void *block_thread(void *arg) {
    os_chan_t *c = arg;
    pthread_mutex_lock(&c->mu);
    while (!c->stop) {
        int64_t b = c->next_b;
        int64_t B = b * c->L - c->V;
        if (B < 0 || (uint64_t)B < c->keep_lo) {     /* gone, or going */
            c->next_b++;
            continue;
        }
        if ((uint64_t)B + (uint64_t)c->N > c->have_hi ||
            b > c->pos_b + OS_CHAN_AHEAD) {
            pthread_cond_wait(&c->cv, &c->mu);       /* not in yet, or far enough ahead */
            continue;
        }
        if (b < c->pos_b - OS_CHAN_BEHIND) {        /* fallen behind: no burst reaches back there */
            c->next_b = c->pos_b - OS_CHAN_BEHIND;
            continue;
        }
        int slot = (int)(b % OS_CHAN_SLOTS);
        if (c->tag[slot] >= b) {
            /* The detector got there first, or has moved on and put a newer
             * block in this slot, which it may be reading: leave it. */
            c->next_b++;
            continue;
        }
        if (c->busy[slot]) {
            pthread_cond_wait(&c->cv, &c->mu);
            continue;
        }
        c->tag[slot] = b;
        c->busy[slot] = 1;
        c->next_b++;
        c->reading = 1;
        c->read_lo = (uint64_t)B;
        pthread_mutex_unlock(&c->mu);

        /* The detector won't overwrite [B, ...) while reading is set
         * (os_chan_overwrite waits), and B >= keep_lo, B + N <= have_hi:
         * all of it is in the ring and stays there for the copy. */
        int complete = c->read(c->read_ctx, (uint64_t)B, (size_t)c->N, c->tfin);

        pthread_mutex_lock(&c->mu);
        c->reading = 0;
        pthread_cond_broadcast(&c->cv);
        pthread_mutex_unlock(&c->mu);

        float complex *dst = &c->slots[(size_t)slot * c->N];
        if (complete)
            fftwf_execute_dft(c->fwd, (fftwf_complex *)c->tfin, (fftwf_complex *)dst);

        pthread_mutex_lock(&c->mu);
        if (complete)
            c->n_fft++;
        else
            c->tag[slot] = -1;
        c->busy[slot] = 0;
        pthread_cond_broadcast(&c->cv);
    }
    pthread_mutex_unlock(&c->mu);
    return NULL;
}

int os_chan_start_thread(os_chan_t *c) {
    if (c->threaded)
        return 0;
    c->tfin = fftwf_alloc_complex(c->N);
    pthread_mutex_init(&c->mu, NULL);
    pthread_cond_init(&c->cv, NULL);
    if (pthread_create(&c->th, NULL, block_thread, c) != 0) {
        pthread_mutex_destroy(&c->mu);
        pthread_cond_destroy(&c->cv);
        fftwf_free(c->tfin);
        c->tfin = NULL;
        return -1;
    }
#ifdef __linux__
    pthread_setname_np(c->th, "blockfft");
#endif
    c->threaded = 1;
    return 0;
}

void os_chan_advance(os_chan_t *c, uint64_t pos, uint64_t have_hi) {
    if (!c->threaded)
        return;
    pthread_mutex_lock(&c->mu);
    c->pos_b = (int64_t)(pos / (uint64_t)c->L);
    c->have_hi = have_hi;
    pthread_cond_broadcast(&c->cv);
    pthread_mutex_unlock(&c->mu);
}

void os_chan_overwrite(os_chan_t *c, uint64_t lo) {
    if (!c->threaded)
        return;
    pthread_mutex_lock(&c->mu);
    if (lo > c->keep_lo)
        c->keep_lo = lo;
    while (c->reading && c->read_lo < lo)
        pthread_cond_wait(&c->cv, &c->mu);
    pthread_mutex_unlock(&c->mu);
}

int os_chan_extract(os_chan_t *c, uint64_t start, uint64_t stop,
                    int bin_shift, float complex *out, int max_out,
                    uint64_t *delay_ns) {
    if (delay_ns)
        *delay_ns = c->delay_ns;
    if (stop <= start || stop - start < (uint64_t)c->ntaps)
        return 0;
    int64_t n_out = (int64_t)(stop - start - c->ntaps + 1) / c->D;
    if (n_out > max_out)
        n_out = max_out;

    const int D = c->D, M = c->M, N = c->N;
    uint64_t m0 = start + c->ntaps - 1;     /* input index under out[0] */
    int tau = (int)(m0 % D);
    const float complex *G = &c->G[(size_t)tau * M];
    int64_t g0 = (int64_t)((m0 - tau) / D);  /* out[0] on the D-sample grid */
    int shift = ((bin_shift % N) + N) % N;

    int64_t k = 0;
    while (k < n_out) {
        int64_t pos = (g0 + k) * D;
        int64_t b = pos / c->L;
        int jj = (int)((pos - b * c->L) / D) + c->V / D;
        int64_t cnt = M - jj;
        if (cnt > n_out - k)
            cnt = n_out - k;

        const float complex *X = get_block(c, b);

        /* Bins shift - M/2 .. shift + M/2 - 1, laid out for the inverse FFT
         * (non-negative k' first), as two runs around the wrap at N. */
        for (int idx = 0; idx < M; idx++) {
            int kp = idx < M / 2 ? idx : idx - M;
            int src = shift + kp;
            if (src < 0) src += N;
            else if (src >= N) src -= N;
            c->zin[idx] = X[src] * G[idx];
        }
        fftwf_execute(c->inv);

        /* Phase reference: the burst's start, not the block's. */
        int64_t B = b * c->L - c->V;
        int64_t rel = (int64_t)start - B;
        int64_t r = ((rel % N) * (int64_t)shift) % N;   /* exact, mod N */
        double ang = 2.0 * M_PI * (double)r / N;
        float complex ph = (float complex)cexp(I * ang);

        for (int64_t i = 0; i < cnt; i++)
            out[k + i] = c->zout[jj + i] * ph;
        k += cnt;
    }
    return (int)n_out;
}
