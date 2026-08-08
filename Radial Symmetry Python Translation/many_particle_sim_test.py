# -*- coding: utf-8 -*-
# many_particle_sim_test.py
"""
Author:   Raghuveer Parthasarathy
Python translation / adaptation of MATLAB many_particle_sim_test.m,
see "TEMP from MATLAB etc/Many particle assessment/"

Description
-----------

Simulates a 2D image containing many randomly-positioned point particles
(via modelNparticleimage.py), finds particles and localizes each using
fo5_rp.py (a partial Python translation of fo5_rp.m: bandpass filtering,
local-max finding, thresholding, and single-particle localization), and
compares each to the known "true" positions to assess localization
accuracy.

Very simple: all particles have the same brightness and size, and there
are no checks for overlapping particles.

Note: the MATLAB version also compared against Gaussian MLE fitting.
fo5_rp.py doesn't implement that option (its MATLAB fitting routine,
gaussfit2DMLE.m, is not part of this repository) -- see fo5_rp.py's
docstring for the full scope of what is and isn't translated. This
script instead compares 'radial', 'centroid', and 'lineargauss'.

Run as a script to reproduce a demo similar to the MATLAB version, or
import assess_accuracy() to use with your own images / fo5_rp() results.

Raghuveer Parthasarathy
See MATLAB many_particle_sim_test.m for revision history
"""

import time

import numpy as np
from scipy.spatial.distance import cdist

from fo5_rp import fo5_rp
from modelNparticleimage import modelNparticleimage


def assess_accuracy(xc_true, yc_true, xc_found, yc_found, label=""):
    """
    Match found particle positions to the nearest true position (by
    Euclidean distance, as in the MATLAB version), and report the mean
    offset and overall (RMS) error, in px.
    """
    Np = len(xc_true)
    Nfound = len(xc_found)
    if Nfound != Np:
        print("*** Warning ***")
        print("  Number of particles found is not equal to the number of simulated")
        print("  particles, probably due to overlapping particle images or a missed")
        print("  detection. Suggest re-running, or adjusting objsize/thresh.")
    if Nfound == 0:
        print(f"{label}\n   No particles found.")
        return np.nan, np.nan, np.nan

    true_pts = np.column_stack((xc_true, yc_true))
    found_pts = np.column_stack((xc_found, yc_found))
    d = cdist(true_pts, found_pts)  # Np x Nfound
    idmap = np.argmin(d, axis=0)  # nearest true particle, for each found one

    dx = xc_found - xc_true[idmap]
    dy = yc_found - yc_true[idmap]
    mx, my = dx.mean(), dy.mean()
    err = np.sqrt(np.mean(dx * dx + dy * dy))

    print(label)
    print(f"   Mean Delta_x = {mx:.3f} px; Mean Delta_y = {my:.3f} px")
    print(f"   Total error (std. of discrepancy) = {err:.3f} px")
    return mx, my, err


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print()
    print("many_particle_sim_test.py:")

    # Parameters -- edit directly, as in the MATLAB script
    N = 9  # each particle image is calculated in an NxN box
    Np = 75  # number of particles
    SNr = 10.0  # signal-to-noise ratio
    Nimage = (600, 800)  # image size (px)
    objsize = (3, 7)  # (bandpass filter size, neighborhood size) -- both
                       # smaller than N, since N was made large to capture
                       # "all" of each simulated particle
    thresh = 0.999  # intensity threshold (percentile option; see fo5_rp.py)

    print(f"Creating model image, {Nimage[0]}x{Nimage[1]} with {Np} particles, "
          f"SNr = {SNr:.1f}...")
    im, xc_true, yc_true = modelNparticleimage(SNr=SNr, Nimage=Nimage, N=N, Np=Np)
    print("  ...done")

    print()
    print("Tracking:")
    print()

    print("Radial-symmetry tracking")
    t0 = time.perf_counter()
    objs_r = fo5_rp(im, "spatialfilter", objsize, thresh, "radial")
    print(f"   ...done.  Elapsed time {time.perf_counter() - t0:.4f} seconds.")

    print("Centroid finding")
    t0 = time.perf_counter()
    objs_c = fo5_rp(im, "spatialfilter", objsize, thresh, "centroid")
    print(f"   ...done.  Elapsed time {time.perf_counter() - t0:.4f} seconds.")

    print("Linear Gaussian (log-intensity parabolic fit) tracking")
    t0 = time.perf_counter()
    objs_lg = fo5_rp(im, "spatialfilter", objsize, thresh, "lineargauss")
    print(f"   ...done.  Elapsed time {time.perf_counter() - t0:.4f} seconds.")

    # Gaussian MLE fitting (as in the MATLAB version) is not implemented
    # here -- see fo5_rp.py's docstring.

    print()
    print("Assessing accuracy:")
    assess_accuracy(xc_true, yc_true, objs_r["x"], objs_r["y"], label="Radial symmetry.")
    assess_accuracy(xc_true, yc_true, objs_c["x"], objs_c["y"], label="Centroid.")
    assess_accuracy(xc_true, yc_true, objs_lg["x"], objs_lg["y"], label="Linear Gaussian.")

    plt.figure(num="simulated image")
    plt.imshow(im, cmap="gray")
    plt.plot(xc_true, yc_true, "go", fillstyle="none", markersize=8, label="true")
    plt.plot(objs_r["x"], objs_r["y"], "r+", label="radial symmetry")
    plt.legend()
    plt.show()
