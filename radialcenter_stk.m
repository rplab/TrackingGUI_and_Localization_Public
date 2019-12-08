% radialcenter_stk.m
%
% Copyright 2011-2012, Raghuveer Parthasarathy, The University of Oregon
%
% Same as radialcenter.m, but allows input of a stack of images
%    Optimized for speed -- array operations to avoid looping through
%    images.
%
% Calculates the center of a 2D intensity distribution.
% Method: Considers lines passing through each half-pixel point with slope
% parallel to the gradient of the intensity at that point.  Considers the
% distance of closest approach between these lines and the coordinate
% origin, and determines (analytically) the origin that minimizes the
% weighted sum of these distances-squared.
% 
% Requires image processing toolbox for use of imfilter; this could easily
% be avoided by a simple 'boxcar' smoothing.
%
%
% Inputs
%   I  : 3D array of 2D intensity distributions (i.e. stack of 2D grayscale images)
%        (Can also be a single 2D image, but this is ~4x slower than 
%         radialcenter.m due to repeated use of repmat)
%
% Outputs
%   xc, yc : the centers of radial symmetry for each 2D slice.
%            xc, yc are each 1 x Nz row vectors, where Nz = #slices
%            px, from px #1 = left/topmost pixel
%            So a shape centered in the middle of a 2*N+1 x 2*N+1
%            square (e.g. from make2Dgaussian.m with x0=y0=0) will return
%            a center value at x0=y0=N+1.
%            (i.e. same convention as gaussfit2Dnonlin.m and my other 
%            particle finding functions)
%            Note that y increases with increasing row number (i.e. "downward")
%   sigma  : Rough measure of the width of the distribution (sqrt. of the 
%            second moment of I - min(I)) for each slice;
%            Not determined by the fit -- output mainly for consistency of
%            formatting compared to my other fitting functions, and to get
%            an estimate of the particle "width."  Can eliminate for speed.
%   meand2 : weighted mean weighted distance^2 from the gradient line distance
%            minimization (Feb. 2013), for each slice.  
%            Not necessary -- output to assess goodness of fit. 
%            Can eliminate for speed.
%
%
% see notes August 19-25, Sept. 9, Sept. 19-20 2011; June 20, 2012
% Raghuveer Parthasarathy
% The University of Oregon
% June 20, 2012 (begun); based on the June 20, 2012 version of
% radialcenter.m
% Nov. 1, 2017.
%    Don't smooth if image is only 3x3
% Sept. 18, 2018
%    Removing weird "mask" -- who put this here?
% Last modified Sept. 18, 2018
%
% Copyright 2011-2017, Raghuveer Parthasarathy
%
%%
% Disclaimer / License  
%   This program is free software: you can redistribute it and/or 
%     modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.
%   This set of programs is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
%   General Public License for more details.
%   You should have received a copy of the GNU General Public License 
%   (gpl.txt) along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

function [xc, yc, sigma, meand2] = radialcenter_stk(I)

% Make sure I is double to allow conv2 etc.
I = double(I);
Nz = size(I,3);  % number of slices

% Number of grid points
Ny = size(I,1); Nx = size(I,2);

% grid coordinates are -n:n, where Nx (or Ny) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5;
xm_onerow = -(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5;
xm = xm_onerow(ones(Ny-1, 1), :);
% similarly replacing
%    ym = repmat((-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)', 1, Nx-1);
ym_onecol = (-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)';  % Note that y increases "downward"
ym = ym_onecol(:,ones(Nx-1,1));


% Calculate derivatives along 45-degree shifted coordinates (u and v)
% Note that y increases "downward" (increasing row number) -- we'll deal
% with this when calculating "m" below.
% Gradient in each slice
dIdu = I(1:Ny-1,2:Nx,:)-I(2:Ny,1:Nx-1,:);
dIdv = I(1:Ny-1,1:Nx-1,:)-I(2:Ny,2:Nx,:);

% Smooth each slice (if the 2D images are > 3px)s
% multidimensional imfilter is faster than repeated calls to conv2; conv2
% is faster for single images
if min([Nx, Ny]>3)
    % Only smooth if image is >3px in the smallest dimension
    h = ones(3)/9;  % simple 3x3 averaging filter
    if Nz>1
        fdu = imfilter(dIdu,h);
        fdv = imfilter(dIdv,h);
    else
        fdu = conv2(dIdu, h, 'same');
        fdv = conv2(dIdv, h, 'same');
    end
else
    fdu = dIdu;
    fdv = dIdv;
end
dImag2 = fdu.*fdu + fdv.*fdv; % gradient magnitude, squared, in each slice

% Slope of the gradient .  Note that we need a 45 degree rotation of
% the u,v components to express the slope in the x-y coordinate system.
% The negative sign "flips" the array to account for y increasing
% "downward"
m = -(fdv + fdu) ./ (fdu-fdv); % gradient slope, in each slice

%    infslope = 9e9;  % replace infinite slope values with this extremely large number
%    m(isinf(m)) = infslope;
infslope = 9e9;  % replace infinite slope values with this extremely large number
m(isinf(m)) = infslope;

% Shorthand "b", which also happens to be the
% y intercept of the line of slope m that goes through each grid midpoint
xmrep = repmat(xm,[1 1 Nz]);
ymrep = repmat(ym,[1 1 Nz]);
b = ymrep - m.*xmrep;

% Weighting: weight by square of gradient magnitude and inverse
% distance to gradient intensity centroid.
sdI2 = squeeze(sum(sum(dImag2,1),2));  % Nz x 1 array of total gradient magnitudes
xcentroid = squeeze(sum(sum(dImag2.*xmrep,1),2))./sdI2;  % Nz x 1 array of x centroids
ycentroid = squeeze(sum(sum(dImag2.*ymrep,1),2))./sdI2;  % Nz x 1 array of y centroids
xcentroidrep = repmat(reshape(xcentroid, [1 1 Nz]), [Ny-1, Nx-1, 1]); 
ycentroidrep = repmat(reshape(ycentroid, [1 1 Nz]), [Ny-1, Nx-1, 1]); 
w = dImag2./sqrt((xmrep-xcentroidrep).*(xmrep-xcentroidrep)+(ymrep-ycentroidrep).*(ymrep-ycentroidrep));

% if the intensity is completely flat, m will be NaN (0/0)
% give these points zero weight (and set m, b = 0 to avoid 0*NaN=NaN)
temp = false(size(m));
temp(isnan(m))=true;
b(temp)=0;
m(temp)=0;
w(temp)=0;

% Calculate the best-fit center for each slice
% m, b, and w above are Ny-1 x Nx-1 x Nz arrays
% Reshape to be two-dimensional ((Ny-1)x(Nx-1)) x Nz arrays, so that each
% column corresponds to one slice
m = reshape(m, [(Ny-1)*(Nx-1), Nz]);
b = reshape(b, [(Ny-1)*(Nx-1), Nz]);
w = reshape(w, [(Ny-1)*(Nx-1), Nz]);

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit_stk (below)
% array operations to calculate each slice's center without looping
[xc, yc] = lsradialcenterfit_stk(m, b, w);

%% calculate mean distance^2 between gradient lines and center, weighted by
% intensity.  A measure of goodness of "fit."  Small == good.
repxc = repmat(xc, size(m,1), 1);
repyc = repmat(yc, size(m,1), 1);
tempd = b-(repyc-m.*repxc);
dmin2 = tempd.*tempd ./ (m.*m+1);  % array of minimal distance-squared values
dImag2 = reshape(dImag2, [(Ny-1)*(Nx-1), Nz]);
meand2 = sum(dmin2.*dImag2,1)./sum(dImag2,1);

%%
% Return output relative to upper left coordinate
xc = xc + (Nx+1)/2.0;
yc = yc + (Ny+1)/2.0;


%% A rough measure of the particle width.
% Not at all connected to center determination, but may be useful for tracking applications; 
% could eliminate for (very slightly) greater speed

% Avoid meshgrid, for greater speed
% [px,py] = meshgrid(1:Nx,1:Ny);
px1row = 1:Nx;    px = px1row(ones(Ny, 1),:);
py1col = (1:Ny)'; py = py1col(:,ones(Nx,1));

% Use intensity at corners as a measure of background
bkgestimate = sum(sum([I(1,1,:) I(1,end,:) I(end,1,:) I(end,end,:)],1),2)/4.0;

Isub = I - repmat(bkgestimate, [Ny, Nx, 1]); %  min(I(:));
Isub(Isub<0) = 0;  % force non-negative; otherwise moment may not make sense
% Subtract minimum from each slice
xoffset = repmat(px,[1 1 Nz]) - repmat(reshape(xc, [1 1 Nz]), [Ny, Nx, 1]);
yoffset = repmat(py,[1 1 Nz]) - repmat(reshape(yc, [1 1 Nz]), [Ny, Nx, 1]);
r2 = xoffset.*xoffset + yoffset.*yoffset;
sigma = squeeze(sqrt(sum(sum(Isub.*r2,1),2)./sum(sum(Isub,1),2))/2); % second moment is 2*Gaussian width
sigma = sigma';  % same orientation as xc, yc


%%

    function [xc, yc] = lsradialcenterfit_stk(m, b, w)
        % least squares solution to determine the radial symmetry center
        % inputs m, b, w are defined on a grid
        % w are the weights for each point
        % each column of the input variables corresponds to an image slice
        
        wm2p1 = w./(m.*m+1);
        sw  = sum(wm2p1,1);
        smmw = sum(m.*m.*wm2p1,1);
        smw  = sum(m.*wm2p1,1);
        smbw = sum(m.*b.*wm2p1,1);
        sbw  = sum(b.*wm2p1,1);
        det = smw.*smw - smmw.*sw;
        xc = (smbw.*sw  - smw.*sbw)./det;    % relative to image center
        yc = (smbw.*smw - smmw.*sbw)./det; % relative to image center
        
    end

end
