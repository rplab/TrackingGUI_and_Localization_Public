% simpleellipsefit.m
% simple function to fit an ellipsoid to a 2D "rod-like" image,
% determining the ellipsoid that matches the covariance matrix of the
% intensity distribution.  Applies a circular 'mask' to avoid bias towards
% the diagonals of the image (optional).
% Determines the ellipse center using simple centroid finding,
% or the user can input the center location (previously found, e.g. via
% radial-symmetry-finding, which works well even for ellipses)
%
% Input: 
%    A : 2D image
%    ctr : (optional) 2-element array with the x and y center locations
%          px, from px #1 = left/topmost pixel
%          If empty, calculate centroid
%    invertopt : (optional)  If true, invert the intensity (max(A(:))-A),
%          before any further calculations. 
%          Default false
%    minsubopt : (optional)  If true, subtract a background value of A 
%          For center determination, the background is simply the minimal
%          value of A.  For ellipse determination, the background is
%          calculated as the mean of A outside the circular mask, or the
%          min. value within the mask (whichever is smaller).
%          default true
%    usemask : (optional)  If true, use a circular mask to avoid bias
%          towards the corners.  Avoids bias, but throws away data
%          (pixels). default true
%    
% Outputs:
%    ctr : 2-element array with the x and y center locations
%          If input, uchanged; else centroid location
%          px, from px #1 = left/topmost pixel
%    theta : orientation, radians; note y increases downward, so "upside
%            down"
%    r : 2-element array with the semimajor and semiminor axis lengths in
%        pixels
%    ecc : eccentricity
%
% Raghuveer Parthasarathy
% August 23, 2012 -- based on previous "quickellipsefit.m"
% Modified Feb. 26, 2013 -- Added usemask option
% Last modified Mar. 14, 2018 -- Changed thresholding to better deal with
% images with high backgrounds

function [ctr, theta, r, ecc] = simpleellipsefit(A, ctr, invertopt, minsubopt, usemask)

% Set default values
if ~exist('ctr', 'var') 
    ctr = [];
end
if ~exist('invertopt', 'var') || isempty(invertopt)
    invertopt = false;
end
if ~exist('minsubopt', 'var') || isempty(minsubopt)
    minsubopt = true;
end
if ~exist('usemask', 'var') || isempty(usemask)
    usemask = true;
end

% Determine the size of the image and generate a grid of coordinates;
% note that y increases downward
A = double(A); % make process-able
[ny, nx] = size(A);
[px, py] = meshgrid(1:nx, 1:ny);

% optional inversion and background subtraction
if invertopt
    % invert intensity
    A = max(A(:)) - A;
end
if minsubopt
    % subtract minimal value, as background
    Asub = A - min(A(:));
else
    Asub = A;
end

% Total intensity (following background subtraction)
sumA = sum(Asub(:));

% Determine center, if not input
if isempty(ctr)
    % centroid
    ctr = [sum(sum(Asub.*px))  sum(sum(Asub.*py))]/sumA;
end

% Determine circular mask in which to examine moments
% consider distance to the center; find the closest approach to the image
% edge.
dx = px-ctr(1);
dy = py-ctr(2);
dr2 = dx.*dx + dy.*dy;
dr2edge = [dr2(1,:) dr2(end,:) dr2(:,1)' dr2(:,end)'];
maskr = true(size(dr2)); 
if usemask
    mindr2edge = min(dr2edge);
    maskr(dr2>mindr2edge) = false;
end

if minsubopt
    % subtract either the value outside the mask or the minimal value 
    % within the mask, as background
    A = A/max(max(A));
    thresh = graythresh(A);
    maskA = A;
    maskA(A<thresh) = 0; 
    % maskA(~maskr) = 0;  % A, background subtracted, with the circular mask applied.  
else
    maskA = A; 
    maskA(~maskr) = 0;  % A, with the circular mask applied.  
end

%figure; imshow(maskA,[]);

% Variance and covariance
newsumA = sum(maskA(:));
xvar = sum(sum(maskA.*dx.*dx))/newsumA;
yvar = sum(sum(maskA.*dy.*dy))/newsumA;
xyvar = sum(sum(maskA.*dx.*dy))/newsumA;

% The covariance matrix is
% c = [xvar xyvar; xyvar yvar];

% Calculate eigenvalues of the variance-covariance matrix
% (These are the *squares* of the major and minor axes of the best-fit ellipse)
D = sqrt((xvar-yvar).*(xvar-yvar) + 4*xyvar*xyvar);
eig1 = 0.5*(xvar+yvar+D);
eig2 = 0.5*(xvar+yvar-D);
r = [sqrt(eig1) sqrt(eig2)];

% Eccentricity
ecc = sqrt(1-(eig2/eig1));  % Note that I'm not squaring the ratio of the 
%   eigenvalues, since these are already the axes^2
% Could also use (eig1-eig2)/(eig1+eig2) as a measure of circularity; I
% think this is the "third eccentricity"

% Angle w.r.t. x-axis.  Note that y increases downward
% Simple calculation of eigenvectors
% Don't need to take square root, since eigenvector angle is the same regardless
% sign1 = sign((eig1-xvar)/xyvar);
% theta = atan(sign1*sqrt(abs((eig1-xvar)/xyvar)));

% look at most sensitive part of atan; borrow code from regionprops
if (yvar > xvar)
    num = yvar - xvar + D;
    den = 2*xyvar;
else
    num = 2*xyvar;
    den = xvar - yvar + D;
end
% Angle w.r.t. x-axis.  Note that y increases downward
theta  = atan(num/den);

% figure; imshow(A,[]);
% hold on
% plot(ctr(1), ctr(2), 'x')
% hold on
% plotphi = linspace(0, 2*pi, 100);
% for j=1:length(plotphi)
%     plot(ctr(1) + 2*r(1)*cos(plotphi)*cos(-theta) + 2*r(2)*sin(plotphi)*sin(-theta),...
%         ctr(2) + 2*r(2)*sin(plotphi)*cos(-theta) - 2*r(1)*cos(plotphi)*sin(-theta), '.');
% end

end

