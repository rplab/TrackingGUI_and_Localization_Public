% radialcenter.m
%
% Copyright 2011-2012, Raghuveer Parthasarathy, The University of Oregon
%
% Calculates the center of a 2D intensity distribution.
% Method: Considers lines passing through each half-pixel point with slope
% parallel to the gradient of the intensity at that point.  Considers the
% distance of closest approach between these lines and the coordinate
% origin, and determines (analytically) the origin that minimizes the
% weighted sum of these distances-squared.
% 
% Does not allow revision / re-weighting of center determination, as this 
% gives negligibly small improvement in accuracy.  See radialcenter_r for
% code.
%
% Inputs
%   I  : 2D intensity distribution (i.e. a grayscale image)
%        Size need not be an odd number of pixels along each dimension
%
% Outputs
%   xc, yc : the center of radial symmetry,
%            px, from px #1 = left/topmost pixel
%            So a shape centered in the middle of a 2*N+1 x 2*N+1
%            square (e.g. from make2Dgaussian.m with x0=y0=0) will return
%            a center value at x0=y0=N+1.
%            (i.e. same convention as gaussfit2Dnonlin.m and my other 
%            particle finding functions)
%            Note that y increases with increasing row number (i.e. "downward")
%   sigma  : Rough measure of the width of the distribution (sqrt. of the 
%            second moment of I - min(I));
%            Not determined by the fit -- output mainly for consistency of
%            formatting compared to my other fitting functions, and to get
%            an estimate of the particle "width."  Can eliminate for speed.
%   meand2 : weighted mean weighted distance^2 from the gradient line distance
%            minimization (Feb. 2013).  
%            Not necessary -- output to assess goodness of fit. 
%            Can eliminate for speed.
%
% see notes August 19-25, Sept. 9, Sept. 19-20 2011
% Raghuveer Parthasarathy
% The University of Oregon
% August 21, 2011 (begun)
% June 15, 2012.
%    Improve speed based on ideas from Eric Corwin -- more efficient
%    checking of infinite slope and undefined slope
% June 20, 2012.  
%    Very minor change: Make "I" double precision (in case it was input as uint16 etc.)
% February 23, 2013.  
%    Very minor change: For sigma calculation, use mean of corners as background
%    Force subtracted image to be non-negative
% Nov. 1, 2017.
%    Don't smooth if image is only 3x3
% Last modified Nov. 1, 2017
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

function [xc, yc, sigma, meand2] = radialcenter(I)

% Make sure I is double to allow conv2 etc.
I = double(I);

% Number of grid points
[Ny, Nx] = size(I);

% grid coordinates are -n:n, where Nx (or Ny) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5;
% The two lines below replace
%    xm = repmat(-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5,Ny-1,1);
% and are faster (by a factor of >15 !)
% -- the idea is taken from the repmat source code
xm_onerow = -(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5;
xm = xm_onerow(ones(Ny-1, 1), :);
% similarly replacing
%    ym = repmat((-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)', 1, Nx-1);
ym_onecol = (-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)';  % Note that y increases "downward"
ym = ym_onecol(:,ones(Nx-1,1));

% Calculate derivatives along 45-degree shifted coordinates (u and v)
% Note that y increases "downward" (increasing row number) -- we'll deal
% with this when calculating "m" below.
dIdu = I(1:Ny-1,2:Nx)-I(2:Ny,1:Nx-1);
dIdv = I(1:Ny-1,1:Nx-1)-I(2:Ny,2:Nx);

% Smoothing -- perhaps should be optional
if min([Nx, Ny]>3)
    % Only smooth if image is >3px in the smallest dimension
    h = ones(3)/9;  % simple 3x3 averaging filter
    fdu = conv2(dIdu, h, 'same');
    fdv = conv2(dIdv, h, 'same');
else
    fdu = dIdu;
    fdv = dIdv;
end
dImag2 = fdu.*fdu + fdv.*fdv; % gradient magnitude, squared

% Slope of the gradient .  Note that we need a 45 degree rotation of 
% the u,v components to express the slope in the x-y coordinate system.
% The negative sign "flips" the array to account for y increasing
% "downward"
m = -(fdv + fdu) ./ (fdu-fdv); 
% m = -(dIdv + dIdu) ./ (dIdu-dIdv); 
% disp('not smoothed!')

infslope = 9e9;  % replace infinite slope values with this extremely large number
m(isinf(m)) = infslope;

% Shorthand "b", which also happens to be the
% y intercept of the line of slope m that goes through each grid midpoint
b = ym - m.*xm;

% Weighting: weight by square of gradient magnitude and inverse 
% distance to gradient intensity centroid.
sdI2 = sum(dImag2(:));
xcentroid = sum(sum(dImag2.*xm))/sdI2;
ycentroid = sum(sum(dImag2.*ym))/sdI2;
w  = dImag2./sqrt((xm-xcentroid).*(xm-xcentroid)+(ym-ycentroid).*(ym-ycentroid));  

% if the intensity is completely flat, m will be NaN (0/0)
% give these points zero weight (and set m, b = 0 to avoid 0*NaN=NaN)
w(isnan(m))=0;
b(isnan(m))=0;
m(isnan(m))=0;


% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit (below)
[xc, yc] = lsradialcenterfit(m, b, w);

%% calculate mean distance^2 between gradient lines and center, weighted by
% intensity.  A measure of goodness of "fit."  Small == good.
tempd = b-(yc-m*xc);
dmin2 = tempd.*tempd ./ (m.*m+1);  % array of minimal distance-squared values
meand2 = sum(sum(dmin2.*dImag2))/sum(dImag2(:));

%%
% Return output relative to upper left coordinate
xc = xc + (Nx+1)/2.0;
yc = yc + (Ny+1)/2.0;

%% A rough measure of the particle width.
% Not at all connected to center determination, but may be useful for tracking applications; 
% could eliminate for (very slightly) greater speed
% Use intensity at corners as a measure of background
bkgestimate = mean([I(1,1) I(1,end) I(end,1) I(end,end)]);
Isub = I - bkgestimate; %  min(I(:));
Isub(Isub<0) = 0;  % force non-negative; otherwise moment may not make sense
% Avoid meshgrid, for greater speed
% [px,py] = meshgrid(1:Nx,1:Ny);
px1row = 1:Nx;    px = px1row(ones(Ny, 1),:);
py1col = (1:Ny)'; py = py1col(:,ones(Nx,1));
 
xoffset = px - xc;
yoffset = py - yc;
r2 = xoffset.*xoffset + yoffset.*yoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width


%% Code for plotting the gradients, for illustration purposes.  Moved to 'illustrate_gradient.m' 
% [h1, h2, h3] = illustrate_gradient(xm, ym, dIdu, dIdv, fdu, fdv, Nx, Ny, xc, yc, m, w, b, I);
% pause

%%

    function [xc, yc] = lsradialcenterfit(m, b, w)
        % least squares solution to determine the radial symmetry center
        
        % inputs m, b, w are defined on a grid
        % w are the weights for each point
        wm2p1 = w./(m.*m+1);
        sw  = sum(sum(wm2p1));
        smmw = sum(sum(m.*m.*wm2p1));
        smw  = sum(sum(m.*wm2p1));
        smbw = sum(sum(m.*b.*wm2p1));
        sbw  = sum(sum(b.*wm2p1));
        det = smw*smw - smmw*sw;
        xc = (smbw*sw - smw*sbw)/det;    % relative to image center
        yc = (smbw*smw - smmw*sbw)/det; % relative to image center
        
    end

end
