% gaussfit2D.m
%
% Two dimensional Gaussian fit
% Uses linear regression to log of function data
% no uncertainty analysis
% fit to form: z = A*exp(-(x-x0)^2 / (2*sigma_x^2))
%                   *exp(-(y-y0)^2 / (2*sigma_y^2))
% * max. value of z must be positive* -- since uses logs to calc. gaussian
% x, y:  1D array.  
% z: 2D array (z at each x, y)
%
% Consider only data > threshold*zmax -- input parameter "threshold" -- if
%    omitted, use threshold = 0.2  Threshold ensures that fit is not
%    dominated by tails of the Gaussian
%
% To use, e.g. x = 1:size(z,1); y = 1:size(z,2); ...
%
% Raghuveer Parthasarathy April 14, 2008

function [A, x0, sigma_x, y0, sigma_y] = gaussfit2D(x, y, z, threshold)

zmax = max(z(:));
if (zmax < 0.0)
    disp('gaussfit.m:  ERROR! zmax must be positive!');
    disp('Press Control-C');
    pause
end
if (nargin < 4)
    threshold = 0.2;
end

Nx = length(x);
Ny = length(y);
% Force x, y to be column, row vectors resp.
if size(x,2)>size(x,1)
    x = x';
end
if size(y,1)>size(y,2)
    y = y';
end

if (length(z(:)) ~= (Nx*Ny))
    disp('Error -- length(z) ~= product of x,y lengths.  Press Control-C');
    pause;
end

% -- Eliminate small or non-positive y-values
x = repmat(x, 1, Ny);
y = repmat(y, Nx, 1);
z = z(:);
x = x(:);
y = y(:);
isgoodz = (z >= threshold*zmax);
z = z(isgoodz);
x = x(isgoodz);
y = y(isgoodz);

logz = log(z);

% MATH
% logz = c1 + c2*x + c3*x^2 + c4*y + c5*y^2, where
% c1 = log(A) - x0^2/(2*sigma_x^2) - y0^2/(2*sigma_y^2)
% c2 = x0 / sigma_x^2
% c3 = -1/(2*sigma_x^2)
% c4 = y0 / sigma_y^2
% c5 = -1/(2*sigma_y^2)


xymat = [ones(length(x),1)  x  x.*x  y  y.*y];
c = xymat \ logz;

sigma_x = sqrt(-1/2/c(3));
sigma_y = sqrt(-1/2/c(5));
x0 = c(2)*sigma_x*sigma_x;
y0 = c(4)*sigma_y*sigma_y;
A = exp(c(1) + x0^2/(2*sigma_x^2) + y0^2/(2*sigma_y^2));


