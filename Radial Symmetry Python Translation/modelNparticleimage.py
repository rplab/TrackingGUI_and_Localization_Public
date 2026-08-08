# -*- coding: utf-8 -*-
# modelNparticleimage.py
"""
Author:   Raghuveer Parthasarathy
Python translation of MATLAB modelNparticleimage.m

Description
-----------

Uses modelImage.py (Simulations of Images and Tracks/modelImage.py) to
create a single image containing several (possibly overlapping)
fluorescent point particles.
Each model particle image is calculated in an NxN box with a random
center position, uniformly distributed over +/- 0.5 px around the box's
central pixel. Each of these NxN arrays is placed in the 2D output array;
pixels not covered by any particle are filled with Poisson background
noise.

Uses 0-indexing of positions (unlike MATLAB)

Inputs
  SNr  : Signal/noise ratio, see modelImage.py. Default 20.0.
  Nimage : output image size (Ny, Nx), px. A single scalar is used for
      both dimensions. Default (200, 200).
  N : size of each particle's CCD image (NxN pixels). Forced to be an
      odd integer to ensure centering in the middle of a pixel.
      Default 7.
  Np : number of particles to simulate. Default 9.
  bkg : background (dark) intensity, Poisson mean. Default 10.0.
  lam : wavelength of light, nm. Default 530.
  NA : numerical aperture. Default 1.3.
  dhr : grid size of the "high resolution" image used for the PSF, nm.
      Default 1.0.

Outputs
  im : simulated 2D image, shape Nimage.
  xc, yc : "true" object centers, px (0-indexed), each length Np.
      y increases with increasing row number (i.e. "downward").

Raghuveer Parthasarathy
See MATLAB modelNparticleimage.m for revision history
"""

import os
import sys

import numpy as np

# modelImage.py lives in the sibling "Simulations of Images and Tracks"
# folder; add it to the path if it isn't already importable.
_SIM_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                         "Simulations of Images and Tracks")
if _SIM_DIR not in sys.path:
    sys.path.insert(0, _SIM_DIR)
from modelImage import modelImage


def modelNparticleimage(SNr=20.0, Nimage=(200, 200), N=7, Np=9, bkg=10.0,
                         lam=530, NA=1.3, dhr=1.0):
    # The main function -- see header comments for details

    if np.isscalar(Nimage):
        Nimage = (Nimage, Nimage)
    Ny, Nx = int(Nimage[0]), int(Nimage[1])

    if N % 2 == 0:
        N = N + 1
        print(f"Forcing N to be odd: {N}")

    # CCD pixel scale; also modelImage.py's default. Only "irrelevant" in
    # that it just sets the units for the random-position generation below.
    scale = 100.0  # nm/px
    maxx0nm = 0.5 * scale  # max sub-pixel displacement, nm

    # Random particle centers (0-indexed px), keeping each NxN box
    # within the image (margin of N px, as in the MATLAB version).
    rng = np.random.default_rng()
    xc = rng.uniform(N, Nx - N, size=Np)
    yc = rng.uniform(N, Ny - N, size=Np)
    xc_round = np.round(xc)
    yc_round = np.round(yc)

    # Sub-pixel offset (nm) of each particle from its box's central pixel;
    # modelImage generates one NxN image per offset in a single call.
    xoffset_nm = (xc - xc_round) * scale
    yoffset_nm = (yc - yc_round) * scale
    im_particles, _, _, _ = modelImage(SNr=SNr, Npx=N, xc=xoffset_nm,
                                        yc=yoffset_nm, bkg=bkg, lam=lam,
                                        NA=NA, scale=scale, dhr=dhr,
                                        maxx0nm=maxx0nm)

    im = np.zeros((Ny, Nx))
    filled = np.zeros((Ny, Nx), dtype=bool)
    half = (N - 1) // 2
    for j in range(Np):
        y0, x0 = int(yc_round[j]), int(xc_round[j])
        yslice = slice(y0 - half, y0 + half + 1)
        xslice = slice(x0 - half, x0 + half + 1)
        im[yslice, xslice] += im_particles[:, :, j]
        filled[yslice, xslice] = True

    # Background: Poisson-distributed noise in pixels not filled in above
    bkgimage = rng.poisson(bkg, size=(Ny, Nx)).astype(float)
    im[~filled] = bkgimage[~filled]

    return im, xc, yc


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    im, xc, yc = modelNparticleimage(SNr=10.0, Nimage=(600, 800), N=9, Np=75)
    plt.figure()
    plt.imshow(im, cmap="gray")
    plt.plot(xc, yc, "r+")
    plt.title("modelNparticleimage.py example")
    plt.show()
