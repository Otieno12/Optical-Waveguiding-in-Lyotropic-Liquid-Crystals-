#!/usr/bin/env python3
"""
Estimate light attenuation along a dye-doped LC tube/fibre from a waveguiding video.

Workflow
  1. Average all frames (suppresses director-fluctuation noise).
  2. Take an intensity profile along the tube axis, averaged across a strip of
     width --strip px, starting at the excitation spot.
  3. Fit  I(x) = A*exp(-k x) + B  (single) and a two-exponential model.
  4. Report 1/e length, attenuation in dB, and lengths in tube diameters.

Usage
  python waveguide_attenuation.py video.avi
  python waveguide_attenuation.py video.avi --scale 0.35 --edges 300 120 340 300
  python waveguide_attenuation.py video.avi --start 469 257 --end 180 140

Options
  --scale   um per pixel (omit to stay in pixels)
  --start   x y of excitation spot (default: brightest smoothed point)
  --end     x y of a point further along the tube (default: PCA of bright region)
  --edges   x1 y1 x2 y2 : two points on opposite tube walls (measures diameter)
  --diameter  tube diameter in the same units as the profile (px if no --scale, else um)
  --fit-from  distance (px) from the spot where the fit starts (default 40)
  --strip   averaging strip width across the tube, px (default 17)
  --length  profile length, px (default 400)
"""
import argparse
import csv
import os

import cv2
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import curve_fit


def mean_frame(path):
    cap = cv2.VideoCapture(path)
    acc, n = None, 0
    while True:
        ok, f = cap.read()
        if not ok:
            break
        g = f[:, :, 0].astype(np.float64) if f.ndim == 3 else f.astype(np.float64)
        acc = g if acc is None else acc + g
        n += 1
    if n == 0:
        raise SystemExit("Could not read any frames from " + path)
    return acc / n, n


def find_axis(img, start=None, end=None):
    sm = gaussian_filter(img, 3)
    if start is None:
        y, x = np.unravel_index(np.argmax(sm), sm.shape)
        start = (float(x), float(y))
    start = np.array(start, float)
    if end is not None:
        d = np.array(end, float) - start
    else:
        mask = gaussian_filter(img, 2) > np.median(img) + 8
        ys, xs = np.nonzero(mask)
        pts = np.c_[xs, ys].astype(float)
        _, _, vt = np.linalg.svd(pts - pts.mean(0), full_matrices=False)
        d = vt[0]
        # point the axis from the spot toward the bulk of the bright region
        if np.dot(pts.mean(0) - start, d) < 0:
            d = -d
    return start, d / np.linalg.norm(d)


def axial_profile(img, start, d, length, strip):
    perp = np.array([-d[1], d[0]])
    L = np.arange(0, length, 1.0)
    s = np.arange(-(strip // 2), strip // 2 + 1)
    prof = []
    for l in L:
        p = start + d * l
        xx = p[0] + perp[0] * s
        yy = p[1] + perp[1] * s
        prof.append(map_coordinates(img, [yy, xx], order=1, mode="nearest").mean())
    return L, np.array(prof)


def single(x, A, k, B):
    return A * np.exp(-k * x) + B


def double(x, A1, k1, A2, k2, B):
    return A1 * np.exp(-k1 * x) + A2 * np.exp(-k2 * x) + B


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("video")
    ap.add_argument("--scale", type=float, default=None)
    ap.add_argument("--start", type=float, nargs=2)
    ap.add_argument("--end", type=float, nargs=2)
    ap.add_argument("--edges", type=float, nargs=4)
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--fit-from", type=float, default=40)
    ap.add_argument("--strip", type=int, default=17)
    ap.add_argument("--length", type=int, default=400)
    ap.add_argument("--out", default=None, help="output prefix")
    a = ap.parse_args()

    img, n = mean_frame(a.video)
    start, d = find_axis(img, a.start, a.end)
    x_px, prof = axial_profile(img, start, d, a.length, a.strip)

    s = a.scale if a.scale else 1.0
    unit = "um" if a.scale else "px"
    x = x_px * s

    # tube diameter
    D = None
    if a.edges:
        x1, y1, x2, y2 = a.edges
        D = np.hypot(x2 - x1, y2 - y1) * s
    elif a.diameter:
        D = a.diameter

    sel = x_px >= a.fit_from
    p1, c1 = curve_fit(single, x[sel], prof[sel], p0=[prof[sel][0], 1 / (100 * s), prof.min()], maxfev=20000)
    k = p1[1]
    e_k = np.sqrt(np.diag(c1))[1]
    res = prof[sel] - single(x[sel], *p1)
    r2 = 1 - np.sum(res ** 2) / np.sum((prof[sel] - prof[sel].mean()) ** 2)

    print("frames averaged      : %d" % n)
    print("spot (x,y) px        : %.0f, %.0f" % tuple(start))
    print("axis direction       : %.3f, %.3f" % tuple(d))
    print("single-exp fit (from %.0f px):" % a.fit_from)
    print("  k                  : %.5f +/- %.5f 1/%s  (R^2 = %.3f)" % (k, e_k, unit, r2))
    print("  1/e length         : %.1f %s" % (1 / k, unit))
    print("  attenuation        : %.4f dB/%s" % (4.343 * k, unit))
    if D:
        print("  tube diameter      : %.1f %s" % (D, unit))
        print("  1/e length / D     : %.2f" % (1 / k / D))
        print("  loss per diameter  : %.2f dB" % (4.343 * k * D))

    try:
        p2, _ = curve_fit(double, x, prof, p0=[prof[0] * 0.6, 0.03 / s, prof[0] * 0.3, 0.005 / s, prof.min()], maxfev=50000)
        fast, slow = sorted([p2[1], p2[3]], reverse=True)
        print("two-exp fit (full profile):")
        print("  fast 1/e = %.1f %s, slow 1/e = %.1f %s" % (1 / fast, unit, 1 / slow, unit))
        print("  slow attenuation   : %.4f dB/%s" % (4.343 * slow, unit))
    except RuntimeError:
        p2 = None
        print("two-exp fit did not converge")

    prefix = a.out or os.path.splitext(a.video)[0]
    with open(prefix + "_profile.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["distance_px", "distance_" + unit, "intensity"])
        for xp, xu, v in zip(x_px, x, prof):
            w.writerow([xp, xu, v])

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].imshow(img, cmap="gray", vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.5))
    end = start + d * a.length
    ax[0].plot([start[0], end[0]], [start[1], end[1]], "c-", lw=1)
    ax[0].plot(*start, "r+", ms=12)
    if a.edges:
        ax[0].plot([a.edges[0], a.edges[2]], [a.edges[1], a.edges[3]], "y-", lw=1.5)
    ax[0].set_title("Mean image + profile line")
    ax[1].plot(x, prof, "k.", ms=3, label="data")
    ax[1].plot(x[sel], single(x[sel], *p1), "r-", label="single exp (1/e = %.0f %s)" % (1 / k, unit))
    if p2 is not None:
        ax[1].plot(x, double(x, *p2), "b--", label="two exp")
    ax[1].set_xlabel("distance (%s)" % unit)
    ax[1].set_ylabel("intensity (a.u.)")
    ax[1].legend()
    ax[2].semilogy(x, np.clip(prof - p1[2], 1e-3, None), "k.", ms=3)
    ax[2].semilogy(x[sel], np.clip(single(x[sel], *p1) - p1[2], 1e-3, None), "r-")
    ax[2].set_xlabel("distance (%s)" % unit)
    ax[2].set_title("background-subtracted (log)")
    plt.tight_layout()
    plt.savefig(prefix + "_attenuation.png", dpi=150)
    print("saved:", prefix + "_profile.csv", prefix + "_attenuation.png")


if __name__ == "__main__":
    main()
