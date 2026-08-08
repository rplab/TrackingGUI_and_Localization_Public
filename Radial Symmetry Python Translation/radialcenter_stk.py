# -*- coding: utf-8 -*-
"""
Particle localization by radial symmetry for image stacks.
Python translation of MATLAB radialcenter_stk.m.

Raghu + Claude Sonnet 5, August 7, 2026


Inputs
------
I : array_like
    3D stack of 2D grayscale images with shape (Ny, Nx, Nz).
    A single 2D image is also accepted and treated as a stack of one.

Outputs
-------
xc, yc : ndarray
    1D arrays of length Nz containing the center of radial symmetry
    for each slice, in pixel coordinates with 0-indexing.
sigma : ndarray
    Rough width estimate for each slice.
meand2 : ndarray
    Weighted mean squared distance for each slice.
"""

import numpy as np


def radialcenter_stk(I):
    I = np.asarray(I, dtype=float)
    if I.ndim == 2:
        I = I[:, :, None]
    elif I.ndim != 3:
        raise ValueError("Input must be a 2D image or a 3D image stack.")

    Ny, Nx, Nz = I.shape

    xm = np.arange(-(Nx - 1) / 2.0 + 0.5, (Nx - 1) / 2.0, 1.0)[None, :]
    ym = np.arange(-(Ny - 1) / 2.0 + 0.5, (Ny - 1) / 2.0, 1.0)[:, None]
    xm = xm[:, : Nx - 1]
    ym = ym[: Ny - 1, :]

    dIdu = I[0:Ny - 1, 1:Nx, :] - I[1:Ny, 0:Nx - 1, :]
    dIdv = I[0:Ny - 1, 0:Nx - 1, :] - I[1:Ny, 1:Nx, :]

    if min(Nx, Ny) > 3:
        pad = ((1, 1), (1, 1), (0, 0))
        dIdu_pad = np.pad(dIdu, pad, mode="constant", constant_values=0)
        dIdv_pad = np.pad(dIdv, pad, mode="constant", constant_values=0)
        fdu = (
            dIdu_pad[0:-2, 0:-2, :]
            + dIdu_pad[0:-2, 1:-1, :]
            + dIdu_pad[0:-2, 2:, :]
            + dIdu_pad[1:-1, 0:-2, :]
            + dIdu_pad[1:-1, 1:-1, :]
            + dIdu_pad[1:-1, 2:, :]
            + dIdu_pad[2:, 0:-2, :]
            + dIdu_pad[2:, 1:-1, :]
            + dIdu_pad[2:, 2:, :]
        ) / 9.0
        fdv = (
            dIdv_pad[0:-2, 0:-2, :]
            + dIdv_pad[0:-2, 1:-1, :]
            + dIdv_pad[0:-2, 2:, :]
            + dIdv_pad[1:-1, 0:-2, :]
            + dIdv_pad[1:-1, 1:-1, :]
            + dIdv_pad[1:-1, 2:, :]
            + dIdv_pad[2:, 0:-2, :]
            + dIdv_pad[2:, 1:-1, :]
            + dIdv_pad[2:, 2:, :]
        ) / 9.0
    else:
        fdu = dIdu
        fdv = dIdv

    dImag2 = fdu * fdu + fdv * fdv
    m = -(fdv + fdu) / (fdu - fdv)
    m[np.isinf(m)] = 9e9

    xm_rep = xm[:, :, None]
    ym_rep = ym[:, :, None]
    b = ym_rep - m * xm_rep

    sdI2 = np.sum(dImag2, axis=(0, 1))
    xcentroid = np.sum(dImag2 * xm_rep, axis=(0, 1)) / sdI2
    ycentroid = np.sum(dImag2 * ym_rep, axis=(0, 1)) / sdI2

    xcentroid = xcentroid[None, None, :]
    ycentroid = ycentroid[None, None, :]
    w = dImag2 / np.sqrt(
        (xm_rep - xcentroid) ** 2 + (ym_rep - ycentroid) ** 2
    )

    nan_mask = np.isnan(m)
    w[nan_mask] = 0
    b[nan_mask] = 0
    m[nan_mask] = 0

    wm2p1 = w / (m * m + 1)
    sw = np.sum(wm2p1, axis=(0, 1))
    smmw = np.sum(m * m * wm2p1, axis=(0, 1))
    smw = np.sum(m * wm2p1, axis=(0, 1))
    smbw = np.sum(m * b * wm2p1, axis=(0, 1))
    sbw = np.sum(b * wm2p1, axis=(0, 1))
    det = smw * smw - smmw * sw
    xc = (smbw * sw - smw * sbw) / det
    yc = (smbw * smw - smmw * sbw) / det

    temp = b - (yc[None, None, :] - m * xc[None, None, :])
    dmin2 = temp * temp / (m * m + 1)
    meand2 = np.sum(dmin2 * dImag2, axis=(0, 1)) / np.sum(dImag2, axis=(0, 1))

    xc = xc + (Nx - 1) / 2.0
    yc = yc + (Ny - 1) / 2.0

    bkgestimate = np.mean(
        np.stack([I[0, 0, :], I[0, -1, :], I[-1, 0, :], I[-1, -1, :]], axis=0),
        axis=0,
    )
    Isub = I - bkgestimate[None, None, :]
    Isub[Isub < 0] = 0

    px = np.arange(Nx)[None, :, None]
    py = np.arange(Ny)[:, None, None]
    xoffset = px - xc[None, None, :]
    yoffset = py - yc[None, None, :]
    r2 = xoffset * xoffset + yoffset * yoffset
    sigma = np.sqrt(np.sum(Isub * r2, axis=(0, 1)) / np.sum(Isub, axis=(0, 1))) / 2.0

    return xc, yc, sigma, meand2
