#!/usr/bin/env python3
"""Functional tests for the hqdn3d VapourSynth plugin.

Usage: python3 test/test_hqdn3d.py [path/to/libhqdn3d.so] [path/to/reference.so]

Without arguments the plugin is expected to be autoloaded (e.g. from the
installed wheel). The optional reference plugin is an older 8 bit only build
whose output the 8 bit path must match bit for bit.

Requires the vapoursynth Python module and numpy.
"""
import sys

import numpy as np
import vapoursynth as vs

core = vs.core

WIDTH, HEIGHT, FRAMES = 96, 64, 12
SEED = 1234


def noise_clip(fmt, depth):
    """A deterministic 8 bit noise pattern scaled to the requested depth."""
    rng = np.random.default_rng(SEED)
    base = core.std.BlankClip(width=WIDTH, height=HEIGHT, length=FRAMES, format=fmt)
    planes = []
    for p in range(base.format.num_planes):
        w = WIDTH >> (base.format.subsampling_w if p else 0)
        h = HEIGHT >> (base.format.subsampling_h if p else 0)
        # smooth gradient + per-frame noise, in 8 bit units
        yy, xx = np.mgrid[0:h, 0:w]
        gradient = (xx * 255 / max(w - 1, 1) + yy * 255 / max(h - 1, 1)) / 2
        frames = []
        for n in range(FRAMES):
            img = gradient + rng.normal(0, 6, size=(h, w)) + 3 * np.sin(n / 2.0)
            frames.append(np.clip(np.rint(img), 0, 255).astype(np.int64) << (depth - 8))
        planes.append(frames)

    dtype = np.uint8 if depth == 8 else np.uint16

    def fill(n, f):
        f = f.copy()
        for p in range(f.format.num_planes):
            np.asarray(f[p])[:] = planes[p][n].astype(dtype)
        return f

    return core.std.ModifyFrame(base, base, fill)


def to_arrays(clip):
    out = []
    for n in range(clip.num_frames):
        f = clip.get_frame(n)
        out.append([np.array(f[p], dtype=np.int64) for p in range(f.format.num_planes)])
    return out


def check(cond, msg):
    if not cond:
        raise SystemExit("FAIL: " + msg)
    print("ok  ", msg)


def main():
    args = sys.argv[1:]
    if args:
        core.std.LoadPlugin(args[0])
    hq = core.hqdn3d.Hqdn3d
    params = dict(lum_spac=4, chrom_spac=3, lum_tmp=6, chrom_tmp=4.5)

    # 1. 8 bit path: bit exact against the reference build, if given. Only
    #    one core exists per process, so the reference runs in a subprocess.
    src8 = noise_clip(vs.YUV420P8, 8)
    out8 = to_arrays(hq(src8, **params))
    if len(args) > 1 and args[1] != "--dump":
        import os
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            dump = os.path.join(tmp, "ref.npz")
            subprocess.check_call([sys.executable, __file__, args[1], "--dump", dump])
            ref8 = np.load(dump)
        same = all(np.array_equal(ref8["f%dp%d" % (n, p)], a)
                   for n, fa in enumerate(out8) for p, a in enumerate(fa))
        check(same, "8 bit output is bit exact with the reference build")
    elif len(args) > 2 and args[1] == "--dump":
        np.savez(args[2], **{"f%dp%d" % (n, p): a for n, fa in enumerate(out8) for p, a in enumerate(fa)})
        return

    # 2. Higher depths filter the same picture to (almost) the same result.
    for fmt, depth, tol in [
        (vs.YUV420P10, 10, 1),
        (vs.YUV420P12, 12, 1),
        (vs.YUV420P16, 16, 1),
    ]:
        out = to_arrays(hq(noise_clip(fmt, depth), **params))
        worst = 0
        for fa, fb in zip(out, out8):
            for a, b in zip(fa, fb):
                down = (a + (1 << (depth - 9))) >> (depth - 8)
                worst = max(worst, int(np.abs(down - b).max()))
        check(worst <= tol, "%d bit output matches 8 bit output within %d (max diff %d)" % (depth, tol, worst))
        check(all(int(a.max()) <= (1 << depth) - 1 for fa in out for a in fa),
              "%d bit output stays within range" % depth)

    # 3. The strength parameters are on the 8 bit scale regardless of the
    #    clip's depth: the same picture must lose the same amount of noise.
    def noise_ratio(fmt, depth):
        src = noise_clip(fmt, depth)
        out = to_arrays(hq(src, lum_spac=6, lum_tmp=8))
        inp = to_arrays(src)
        hdiff = lambda a: np.std(np.diff(a.astype(np.float64), axis=1))
        return hdiff(out[-1][0]) / hdiff(inp[-1][0])

    r8 = noise_ratio(vs.GRAY8, 8)
    r16 = noise_ratio(vs.GRAY16, 16)
    check(r8 < 0.9, "8 bit clip is denoised (residual noise %.3f)" % r8)
    check(abs(r8 - r16) < 0.01, "16 bit clip is denoised as strongly as 8 bit (%.3f vs %.3f)" % (r16, r8))

    # 3b. From 12 bit on the coefficient tables resolve one source LSB, so
    #     noise far below an 8 bit step is still filtered at 16 bit.
    rng = np.random.default_rng(SEED)
    base = core.std.BlankClip(width=WIDTH, height=HEIGHT, length=6, format=vs.GRAY16, color=[32768])
    fine = [np.clip(32768 + rng.normal(0, 4, (HEIGHT, WIDTH)), 0, 65535).astype(np.uint16) for _ in range(6)]

    def fill_fine(n, f):
        f = f.copy()
        np.asarray(f[0])[:] = fine[n]
        return f

    fine_src = core.std.ModifyFrame(base, base, fill_fine)
    fine_out = to_arrays(hq(fine_src, lum_spac=4, lum_tmp=6))
    fine_std = float(fine_out[-1][0].astype(np.float64).std())
    check(fine_std < 1.0, "16 bit sub-LSB noise is filtered (std 4.0 -> %.2f)" % fine_std)

    # 4. Format coverage and constant clip pass-through.
    for fmt in [vs.GRAY8, vs.GRAY10, vs.GRAY16, vs.YUV444P10, vs.YUV422P12, vs.YUV420P14, vs.YUV444P16]:
        f = core.get_video_format(fmt)
        value = 1 << (f.bits_per_sample - 1)
        clip = core.std.BlankClip(width=WIDTH, height=HEIGHT, length=4, format=fmt, color=[value] * f.num_planes)
        out = to_arrays(hq(clip, **params))
        check(all(np.all(a == value) for fa in out for a in fa), "%s constant clip is unchanged" % f.name)

    # 5. Out of range samples must not crash (they index the LUT).
    hot = core.std.BlankClip(width=WIDTH, height=HEIGHT, length=3, format=vs.GRAY10, color=[1023])

    def overflow(n, f):
        f = f.copy()
        np.asarray(f[0])[:] = 4000
        return f

    hot = core.std.ModifyFrame(hot, hot, overflow)
    out = to_arrays(hq(hot, **params))
    check(all(int(a.max()) == 1023 for fa in out for a in fa), "out of range 10 bit samples are clamped")

    # 6. Seeking: requesting frames out of order works at 16 bit.
    clip = hq(noise_clip(vs.YUV420P16, 16), **params)
    for n in [7, 2, 11, 0, 5]:
        clip.get_frame(n)
    check(True, "out of order frame requests at 16 bit")

    # 7. Unsupported formats are rejected.
    for fmt in [vs.RGB24, vs.GRAYS, vs.YUV444PS]:
        try:
            hq(core.std.BlankClip(format=fmt))
        except vs.Error:
            check(True, "%s is rejected" % core.get_video_format(fmt).name)
        else:
            raise SystemExit("FAIL: %s was accepted" % core.get_video_format(fmt).name)

    print("all tests passed")


if __name__ == "__main__":
    main()
