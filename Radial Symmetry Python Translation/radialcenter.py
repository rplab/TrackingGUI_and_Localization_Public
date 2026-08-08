# -*- coding: utf-8 -*-
# radialcenter.py
"""
Author:   Raghuveer Parthasarathy
Created on Mon Oct 31, 2022
Last modified August 7, 2026 (Raghu + Claude Sonnet 5 -- update meshgrid)

Description
-----------

Particle localization by radial symmetry
Python translation of MATLAB radialcenter.m
Uses 0 indexing of positions (unlike MATLAB)
NOTE: *Does not* optimize for image stacks (like radialcenter_stk.m);
      just single image

Copyright 2011-2022, Raghuveer Parthasarathy, The University of Oregon

Calculates the center of a 2D intensity distribution.
Method: Considers lines passing through each half-pixel point with slope
parallel to the gradient of the intensity at that point.  Considers the
distance of closest approach between these lines and the coordinate
origin, and determines (analytically) the origin that minimizes the
weighted sum of these distances-squared.
Applies simple smoothing if size > 3x3

Inputs
  I  : 2D intensity distribution (i.e. a grayscale image)
       Size need not be an odd number of pixels along each dimension

Outputs
  xc, yc : the center of radial symmetry,
           px, from px #0 = left/topmost pixel
           So a shape centered in the middle of a 2*N+1 x 2*N+1
           square will return a center value at x0=y0=N. (Unlike MATLAB: N+1)
           Note that y increases with increasing row number (i.e. "downward")
  sigma  : Rough measure of the width of the distribution (sqrt. of the 
           second moment of I - min(I));
           Not determined by the fit -- output mainly for consistency of
           formatting compared to my other fitting functions, and to get
           an estimate of the particle "width."  Can eliminate for speed.
  meand2 : weighted mean weighted distance^2 from the gradient line distance
           minimization (Feb. 2013).  
           Not necessary -- output to assess goodness of fit. 
           Can eliminate for speed.

To do:
    - Test more (like MATLAB version)

see notes August 19-25, Sept. 9, Sept. 19-20 2011
Raghuveer Parthasarathy
The University of Oregon
August 21, 2011 (begun)
   Very minor change: For sigma calculation, use mean of corners as background
   Force subtracted image to be non-negative
%%
Disclaimer / License  
  This program is free software: you can redistribute it and/or 
    modify it under the terms of the GNU General Public License as 
    published by the Free Software Foundation, either version 3 of the 
    License, or (at your option) any later version.
  This set of programs is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
  General Public License for more details.
  You should have received a copy of the GNU General Public License 
  (gpl.txt) along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%


"""

import numpy as np   # Will assume numpy is already imported as np !


def radialcenter(I):
    # The main function -- see header comments for details
    
    I = np.asarray(I, dtype=float)
    (Ny, Nx) = I.shape
    # grid coordinates are -n:n, where Nx (or Ny) = 2*n+1
    # grid midpoint coordinates are -n+0.5:n-0.5
    # Use 1D coordinate vectors with broadcasting instead of meshgrid,
    # similar to the faster MATLAB repmat approach.
    xm = np.arange(-(Nx-1)/2.0 + 0.5, (Nx-1)/2.0, 1.0)[None, :]
    ym = np.arange(-(Ny-1)/2.0 + 0.5, (Ny-1)/2.0, 1.0)[:, None]
    # Calculate derivatives along 45-degree shifted coordinates (u and v)
    # Note that y increases "downward" (increasing row number) -- we'll deal
    # with this when calculating "m" below.
    dIdu = I[0:Ny-1,1:]-I[1:,0:Nx-1]
    dIdv = I[0:Ny-1,0:Nx-1]-I[1:,1:]
    # Smoothing -- perhaps should be optional
    fdu = dIdu # will overwrite if smoothing
    fdv = dIdv
    if np.min((Nx, Ny))>3:
        # Only smooth if image is >3px in the smallest dimension
        # Smooth by simple 3x3 boxcar, which I'll code directly rather than
        #    calling a convolution.
        # Zero-pad (expand by 1 on each side), then sum the 9 shifted views
        #    of the padded array -- equivalent to a 3x3 boxcar but vectorized
        #    (avoids a slow pixel-by-pixel Python loop)
        dIdu_pad = np.pad(dIdu, 1) # dIdu array is size Ny-1, Nx-1
        dIdv_pad = np.pad(dIdv, 1) # dIdv array is size Ny-1, Nx-1
        fdu = sum(dIdu_pad[j:j+Ny-1, k:k+Nx-1] for j in range(3) for k in range(3)) / 9.0
        fdv = sum(dIdv_pad[j:j+Ny-1, k:k+Nx-1] for j in range(3) for k in range(3)) / 9.0
    dImag2 = fdu*fdu + fdv*fdv # gradient magnitude, squared
    
    # Slope of the gradient .  Note that we need a 45 degree rotation of 
    # the u,v components to express the slope in the x-y coordinate system.
    # The negative sign "flips" the array to account for y increasing
    # "downward"
    m = -(fdv + fdu) / (fdu-fdv)
    # Not smoothed version: m = -(dIdv + dIdu) ./ (dIdu-dIdv)
    
    infslope = 9e9 #replace infinite slope values with this extremely large number
    m[np.isinf(m)] = infslope
    
    # Shorthand "b", which also happens to be the
    # y intercept of the line of slope m that goes through each grid midpoint
    b = ym - m*xm
    
    # Weighting: weight by square of gradient magnitude and inverse 
    # distance to gradient intensity centroid.
    sdI2 = np.sum(dImag2)
    xcentroid = np.sum(dImag2*xm)/sdI2
    ycentroid = np.sum(dImag2*ym)/sdI2
    w  = dImag2/np.sqrt((xm-xcentroid)*(xm-xcentroid) + 
                        (ym-ycentroid)*(ym-ycentroid))
    
    # if the intensity is completely flat, m will be NaN (0/0)
    # give these points zero weight (and set m, b = 0 to avoid 0*NaN=NaN)
    w[np.isnan(m)]=0
    b[np.isnan(m)]=0
    m[np.isnan(m)]=0
    
    # least-squares minimization to determine the translated coordinate
    # system origin (xc, yc) such that lines y = mx+b have
    # the minimal total distance^2 to the origin:
    # Unilke the MATLAB version, where I have a separate function
    #   for this (lsradialcenterfit), I'll just write the calculation here:
    # Note m, b, w are defined on a grid;  w are the weights for each point
    wm2p1 = w/(m*m+1)
    sw  = np.sum(wm2p1)
    smmw = np.sum(m*m*wm2p1)
    smw  = np.sum(m*wm2p1)
    smbw = np.sum(m*b*wm2p1)
    sbw  = np.sum(b*wm2p1)
    det = smw*smw - smmw*sw
    xc = (smbw*sw - smw*sbw)/det # relative to image center
    yc = (smbw*smw - smmw*sbw)/det # relative to image center
    
    # calculate mean distance^2 between gradient lines and center, weighted by
    # intensity.  A measure of goodness of "fit."  Small == good.
    tempd = b-(yc-m*xc)
    dmin2 = tempd*tempd / (m*m+1)  # array of minimal distance-squared values
    meand2 = np.sum(dmin2*dImag2)/np.sum(dImag2)
    
    # Return output relative to upper left coordinate
    xc = xc + (Nx-1)/2.0
    yc = yc + (Ny-1)/2.0
    
    # A rough measure of the particle width.
    # Not at all connected to center determination, 
    #   but may be useful for tracking applications.
    # could eliminate for (very slightly) greater speed
    # Use intensity at corners as a measure of background
    bkgestimate = np.mean((I[0,0], I[0,-1], I[-1,0], I[-1,-1]))
    Isub = I - bkgestimate
    Isub[Isub<0] = 0 # force non-negative; otherwise moment may not make sense
    # Avoid meshgrid to reduce temporary arrays
    px = np.arange(Nx)[None, :]
    py = np.arange(Ny)[:, None]
    xoffset = px - xc
    yoffset = py - yc
    r2 = xoffset*xoffset + yoffset*yoffset
    
    # second moment is 2*Gaussian width
    sigma = np.sqrt(np.sum(Isub*r2)/np.sum(Isub))/2
    
    return xc, yc, sigma, meand2

