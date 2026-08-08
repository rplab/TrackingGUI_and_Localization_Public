# -*- coding: utf-8 -*-
# fo5_rp.py
"""
Author:   Raghuveer Parthasarathy
Partial Python translation of MATLAB fo5_rp.m
(see "TEMP from MATLAB etc/fo5_rp.m")

Description
-----------

Finds objects (particles) in a 2D image: determines neighborhoods likely
to contain one particle each (by bandpass-filtering, finding local
maxima, and thresholding), then refines each position with a
single-particle localization method.

Uses 0-indexing of positions (unlike MATLAB).

Scope of this translation
--------------------------
fo5_rp.m supports several processing/fitting/orientation options that
call other MATLAB functions (calcthreshpts.m, gradientvote.m,
gaussfit2DMLE.m, gaussfit2Dnonlin.m, gaussfit2D.m, gauss2dcirc.m,
simpleellipsefit.m, and others) that are not part of this repository.
This translation implements only the options with available Python
(or directly-portable) implementations:
  processopt  : 'spatialfilter' (via bpass.py) or 'none'.
                'gradientvote' raises NotImplementedError.
  thresh      : all three options from the MATLAB docstring are
                implemented (percentile, std.-dev.-above-median, and
                brightest-N).
  fitstr      : 'radial' (radialcenter.py / radialcenter_stk.py),
                'centroid' (unweighted intensity centroid, no background
                subtraction, as in the MATLAB version), and
                'lineargauss' (parabolic fit to log-intensity; adapted
                from fo4_rp.m's inline implementation, since fo5_rp.m's
                own gaussfit2D.m is not available here).
                Other fitstr values raise NotImplementedError.
  Orientation-finding (ellipse fitting) is not implemented.
  nhoodctrs   : supported (skips local-max finding).

Inputs
  img : image (2D array) to locate objects within
  processopt : 'spatialfilter' [default] or 'none'
  processparam : for 'spatialfilter', (bpfiltsize, nsize) -- or a single
      value used for both. "0" bpfiltsize means no filtering.
      For 'none', nsize (the neighborhood size).
      nsize is used to build a disk structuring element (matching
      MATLAB's strel('disk', floor(nsize/2), 0)) for local-max finding
      and for the single-particle neighborhood/crop size.
  thresh : intensity threshold. See fo5_rp.m header for the three
      options (percentile / std-dev / brightest-N), reproduced above.
  fitstr : 'radial' [default], 'centroid', or 'lineargauss'
  try1pernhood : unused placeholder (this translation's local-max
      finding is always effectively "one per neighborhood", via the
      disk-shaped dilation footprint) -- kept for interface parity.
  nhoodctrs : optional (N, 2) array of (x, y) neighborhood centers, px,
      0-indexed. If given, skips local-max finding.

Output
  A dict of 1D arrays (one entry per found object), 0-indexed px:
      'x', 'y' : position
      'mass'   : background-subtracted integrated intensity within the
                 disk neighborhood
      'id'     : 1, 2, 3, ...
      'sigma'  : width estimate (see fitstr above; 0 for 'centroid')
      'meand2' : goodness-of-fit (radial only; 0 otherwise)
  (MATLAB's frame/trackid rows are omitted -- set those by the caller.)

Raghuveer Parthasarathy
See MATLAB fo5_rp.m for revision history
"""

import numpy as np
from scipy.ndimage import maximum_filter

from bpass import bpass
from radialcenter import radialcenter
from radialcenter_stk import radialcenter_stk


def _disk_footprint(radius):
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (xx * xx + yy * yy) <= radius * radius


def _centroid_fit(crop):
    # Intensity centroid; no background subtraction (as in fo5_rp.m)
    Ny, Nx = crop.shape
    y, x = np.mgrid[0:Ny, 0:Nx]
    total = np.sum(crop)
    xc = np.sum(crop * x) / total
    yc = np.sum(crop * y) / total
    return xc, yc


def _lineargauss_fit(crop):
    # Parabolic (linear-regression) fit to log(intensity):
    #   log(I) = a0 + a1*x^2 + a2*x + a3*y^2 + a4*y
    # Adapted from fo4_rp.m's inline 'gauss' fit (fo5_rp.m's own
    # gaussfit2D.m helper is not available in this repository).
    Ny, Nx = crop.shape
    y, x = np.mgrid[0:Ny, 0:Nx]
    positive = crop > 0
    xa = x[positive].astype(float)
    ya = y[positive].astype(float)
    logint = np.log(crop[positive])
    A = np.column_stack((np.ones_like(xa), xa * xa, xa, ya * ya, ya))
    coeffs, *_ = np.linalg.lstsq(A, logint, rcond=None)
    _, a1, a2, a3, a4 = coeffs
    xcent = -a2 / (2.0 * a1) if a1 != 0 else np.nan
    ycent = -a4 / (2.0 * a3) if a3 != 0 else np.nan
    sigma_x2 = -1.0 / (2.0 * a1) if a1 < 0 else 0.0
    sigma_y2 = -1.0 / (2.0 * a3) if a3 < 0 else 0.0
    sigma = 0.5 * (np.sqrt(sigma_x2) + np.sqrt(sigma_y2))
    return xcent, ycent, sigma


def fo5_rp(img, processopt="spatialfilter", processparam=7, thresh=0.999,
           fitstr="radial", try1pernhood=False, nhoodctrs=None):
    # The main function -- see header comments for details

    img = np.asarray(img, dtype=float)
    Ny, Nx = img.shape

    processopt_l = (processopt or "none").lower()
    processparam_arr = np.atleast_1d(processparam)
    if processopt_l == "spatialfilter":
        bpfiltsize = float(processparam_arr[0])
        nsize = float(processparam_arr[-1])
        if bpfiltsize < 1:
            processedimg = img
        else:
            processedimg = bpass(img, 1.0, bpfiltsize)
    elif processopt_l == "none":
        nsize = float(processparam_arr[0])
        processedimg = img
    elif processopt_l == "gradientvote":
        raise NotImplementedError(
            'processopt="gradientvote" requires gradientvote.m, which is '
            "not available in this repository.")
    else:
        raise ValueError(f'Unknown processopt "{processopt}"')

    fitstr_l = fitstr.lower()
    if fitstr_l not in ("radial", "centroid", "lineargauss"):
        raise NotImplementedError(
            f'fitstr="{fitstr}" is not implemented in this Python '
            "translation (its MATLAB fitting routine is not part of this "
            'repository). Supported: "radial", "centroid", "lineargauss".')

    # Threshold option -- see header comments
    if thresh >= 1.0:
        threshopt = 3
    elif thresh >= 0.0:
        threshopt = 1
    else:
        threshopt = 2
        thresh = -thresh

    radius = int(np.floor(nsize / 2.0))
    lenx = 2 * radius + 1
    footprint = _disk_footprint(radius)

    if nhoodctrs is None:
        dilated = maximum_filter(processedimg, footprint=footprint, mode="nearest")
        is_local_max = processedimg == dilated
        if threshopt == 1:
            cutoff = np.percentile(processedimg, thresh * 100.0)
            rows, cols = np.nonzero(is_local_max & (processedimg > cutoff))
        elif threshopt == 2:
            # ddof=1 to match MATLAB's std() (sample, not population, std)
            cutoff = np.median(processedimg) + thresh * np.std(processedimg, ddof=1)
            rows, cols = np.nonzero(is_local_max & (processedimg > cutoff))
        else:  # threshopt == 3: brightest N candidates
            rows, cols = np.nonzero(is_local_max)
            n_keep = int(round(thresh))
            if len(rows) > n_keep:
                order = np.argsort(processedimg[rows, cols])[::-1][:n_keep]
                rows, cols = rows[order], cols[order]
        x = cols.astype(float)
        y = rows.astype(float)
    else:
        nhoodctrs = np.asarray(nhoodctrs, dtype=float)
        x = nhoodctrs[:, 0]
        y = nhoodctrs[:, 1]
        good = ~(np.isnan(x) | np.isnan(y))
        x, y = x[good], y[good]

    # Discard candidates too close to the edge for a full neighborhood
    x_round = np.round(x).astype(int)
    y_round = np.round(y).astype(int)
    keep = ((x_round - radius >= 0) & (x_round + radius < Nx) &
            (y_round - radius >= 0) & (y_round + radius < Ny))
    x_round, y_round = x_round[keep], y_round[keep]

    n = len(x_round)
    if n == 0:
        empty = np.array([])
        return {"x": empty, "y": empty, "mass": empty, "id": empty,
                "sigma": empty, "meand2": empty}

    # Crop neighborhoods from the *original* (unfiltered) image
    crops = np.empty((lenx, lenx, n))
    for k in range(n):
        r0, c0 = y_round[k] - radius, x_round[k] - radius
        crops[:, :, k] = img[r0:r0 + lenx, c0:c0 + lenx]

    # Brightness ("mass"): sum within the disk neighborhood, background-
    # subtracted using the mean of the pixels outside the disk
    mass = np.empty(n)
    for k in range(n):
        crop = crops[:, :, k]
        bkg = crop[~footprint].mean()
        mass[k] = crop[footprint].sum() - footprint.sum() * bkg

    sigma = np.zeros(n)
    meand2 = np.zeros(n)
    if fitstr_l == "radial":
        if n < 10:
            xcent = np.empty(n)
            ycent = np.empty(n)
            for k in range(n):
                xcent[k], ycent[k], sigma[k], meand2[k] = radialcenter(crops[:, :, k])
        else:
            xcent, ycent, sigma, meand2 = radialcenter_stk(crops)
        # Sanity check (as in fo5_rp.m): replace wild fits with the centroid
        center = (lenx - 1) / 2.0
        bad = (np.abs(xcent - center) > 1.5 * lenx) | (np.abs(ycent - center) > 1.5 * lenx)
        for k in np.nonzero(bad)[0]:
            xcent[k], ycent[k] = _centroid_fit(crops[:, :, k])
    elif fitstr_l == "centroid":
        xcent = np.empty(n)
        ycent = np.empty(n)
        for k in range(n):
            xcent[k], ycent[k] = _centroid_fit(crops[:, :, k])
        # sigma left at 0, per fo5_rp.m's documented contract for 'centroid'
    else:  # lineargauss
        xcent = np.empty(n)
        ycent = np.empty(n)
        for k in range(n):
            xcent[k], ycent[k], sigma[k] = _lineargauss_fit(crops[:, :, k])

    xn = xcent + (x_round - radius)
    yn = ycent + (y_round - radius)

    return {"x": xn, "y": yn, "mass": mass, "id": np.arange(1, n + 1),
            "sigma": sigma, "meand2": meand2}
