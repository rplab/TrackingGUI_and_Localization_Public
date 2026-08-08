# -*- coding: utf-8 -*-
# bpass.py
"""
Author:   Raghuveer Parthasarathy
Python translation of bpass.m (J. Crocker / D. Grier / E. Weeks;
long a standard part of the Grier-lab particle-tracking toolchain, and
called by fo4_rp.m / fo5_rp.m for image pre-processing before finding
particle candidates).

Description
-----------

Implements a real-space bandpass filter that suppresses pixel noise and
long-wavelength image variations while retaining information of a
characteristic size ("unsharp masking"): convolve with a small Gaussian
to suppress noise, convolve with a boxcar to estimate the local
background, and subtract.

Inputs
  image_array : 2D array to filter
  lnoise : characteristic lengthscale of noise, px. Additive noise
      averaged over this length should vanish. May be 0, in which case
      only the boxcar background subtraction is performed.
  lobject : integer length, px, somewhat larger than a typical object.
      May be 0 (or falsy), in which case only the Gaussian low-pass
      "blurring" is done, without background subtraction.
  threshold : after filtering, pixels below this value are set to 0.
      Default 0.0.

Output
  res : filtered image, same shape as image_array

Raghuveer Parthasarathy
Python translation of bpass.m
"""

import numpy as np
from scipy.ndimage import convolve1d


def bpass(image_array, lnoise, lobject, threshold=0.0):
    image_array = np.asarray(image_array, dtype=float)

    if lnoise == 0:
        gaussian_kernel = np.array([1.0])
    else:
        half_g = int(np.ceil(5 * lnoise))
        xs = np.arange(-half_g, half_g + 1)
        gaussian_kernel = np.exp(-(xs / (2.0 * lnoise)) ** 2)
        gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()

    use_boxcar = bool(lobject)
    if use_boxcar:
        half_b = int(round(lobject))
        boxcar_kernel = np.ones(2 * half_b + 1)
        boxcar_kernel = boxcar_kernel / boxcar_kernel.sum()

    def _conv2sep(img, kernel):
        # Separable 2D convolution (kernel is symmetric, so correlation ==
        # convolution); zero-padded, matching MATLAB's conv2(...,'same').
        out = convolve1d(img, kernel, axis=1, mode="constant", cval=0.0)
        out = convolve1d(out, kernel, axis=0, mode="constant", cval=0.0)
        return out

    gconv = _conv2sep(image_array, gaussian_kernel)
    if use_boxcar:
        bconv = _conv2sep(image_array, boxcar_kernel)
        filtered = gconv - bconv
    else:
        filtered = gconv

    # Zero out the edges (not enough support for a full-width convolution)
    lzero = int(round(max(lobject, np.ceil(5 * lnoise))))
    if lzero > 0:
        filtered[:lzero, :] = 0
        filtered[-lzero:, :] = 0
        filtered[:, :lzero] = 0
        filtered[:, -lzero:] = 0

    filtered[filtered < threshold] = 0
    return filtered
