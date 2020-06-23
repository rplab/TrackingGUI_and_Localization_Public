% gaussfit2DMLE.m
%
% function to fit a 2D (symmetric) Gaussian using Maximum-Likelihood
% Estimation (MLE).
% Fits to a Gaussian with Poisson-distributed noise and a constant (zero
% standard deviation) background
% MLE equations from Abraham et al., Opt. Express, 17:23352 (2009).
% checks if fit center is outside of image; if so, output x0=y0=center.
% (MLE can give nonsensical results, esp. at low signal/noise)
%
% Input: 
%    z : 2D array (e.g. particle image). Can be a 1D array if x and y pixel
%        values are separately input (px, py, below)
%    tolerance in z for fitting (default 1e-6 if empty)  
%    params0 : [optional] starting values of the parameters
%              If not input or empty, default values
%              1 - constant offset (minimal value of z)
%              2 - x-coordinate of the "center of mass" of the 4 brightest pixels
%              3 - y-coordinate of the "center of mass" of the 4 brightest pixels
%              4 - sigma   (quarter of image, roughly)
%              5 - amplitude (max. z - min. z)
%    px, py : x and y pixel values (optional -- usually omit!).  Useful if,
%        e.g., fitting just part of an image.  If these are input, params0
%        must be input, since default calculation won't work.
% Output:
%    A  : Gaussian amplitude
%    x0, y0 : Gaussian center (in px, from px #1 = left/topmost pixel)
%             So a Gaussian centered in the middle of a 2*N+1 x 2*N+1
%             square (e.g. from make2Dgaussian.m with x0=y0=0) will return
%             a center value at x0=y0=N+1.
%             Note that y increases "downward" (with increasing row no.)
%    sigma: Std dev. of Gaussian (assumed same for x, y)
%    offset : constant offset
%
% Raghuveer Parthasarathy 
% February 7, 2012
% March 15, 2012 -- px, py inputs
% April 6, 2012 -- check if center is outside the image
% Sept. 19, 2017: minor fix to offset calculation, out-of-bounds calculation.
% last modified June 11, 2020: minor fix to out-of-bounds calculation 

function [A, x0, y0, sigma, offset] = gaussfit2DMLE(z, tolz, params0, px, py)

if ~exist('px', 'var') || isempty(px)
    [ny,nx] = size(z);
    [px,py] = meshgrid(1:nx,1:ny);
else
    if ~exist('params0', 'var')
        errordlg('Error in gaussfit2DMLE: params0 must be input if px, py are input')
        % There will be errors in determining default starting parameter
        % values.  Note that ny, nx will be undefined
    end
end

% defaults for initial parameter values, and lower and upperbounds
if ~exist('tolz', 'var') || isempty(tolz)
    tolz = 1e-6;
end
if ~exist('params0', 'var') || isempty(params0)
    % Default: "center of mass" of the 3x3 square centered on the brightest pixel 
    [maxzy, maxzx] = find(z==max(z(:)));
    % If there are multiple equal-brightness maxima (e.g. if the image is
    % saturated) just use the image center:
    if size(maxzx,1)>1
        maxzx = ceil(nx/2);
        maxzy = ceil(ny/2);
    end
    % shift if the brightest pixel is at the edge
    if maxzx < 2, maxzx = 2; end
    if maxzx > (nx-1), maxzx = nx-1; end
    if maxzy <2, maxzy = 2; end
    if maxzy > (ny-1), maxzy = ny-1; end
    zc = z(maxzy-1:maxzy+1,maxzx-1:maxzx+1);
    % if zc isn't flat, calculate "center of mass"
    if max(abs(diff(zc(:))))>(100*eps)
        zc = zc-min(zc(:));  % subtract the minimal value as an "offset"
        pxc = px(maxzy-1:maxzy+1,maxzx-1:maxzx+1);
        pyc = py(maxzy-1:maxzy+1,maxzx-1:maxzx+1);
        sumzc = sum(zc(:));
        params0_2 = sum(sum(zc.*pxc))/sumzc;
        params0_3 = sum(sum(zc.*pyc))/sumzc;
    else
        params0_2 = maxzx;  % simple midpoint
        params0_3 = maxzy;  % simple midpoint
    end
    params0 = [min(z(:)), params0_2, params0_3, min([nx ny])/4, max(z(:))-min(z(:))];
end

% More fitting options
fminoptions.TolFun = tolz;  %  % MATLAB default is 1e-6
fminoptions.TolX = 1e-6';  % default is 1e-6
fminoptions.Display = 'off'; % 'off' or 'final'; 'iter' for display at each iteration
fminoptions.LargeScale = 'off';  % use the medium-scale algorithm,
   % since I'm not supplying the gradient
   
% Add an offset, to avoid the chance of small numbers causing problems for
% the estimator function (only makes a difference for small signal/noise,
% less than about 6.)
zoff = 10.0;  % Does choice of zoff matter?
z = z+zoff;
params0(1) = params0(1) + zoff;

params = fminunc(@(P) objfun(P,px,py,z),params0,fminoptions);
A = params(5);
x0 = params(2);
y0 = params(3);
sigma = params(4);
offset = params(1) - zoff;

% Check if particle center is outside of image; return simple centroid if so
if (x0<min(px(:))) || (x0>max(px(:))) || (x0<min(py(:))) || (x0>max(py(:)))
    sz = sum(z(:));
    x0 = sum(sum(z.*px))/sz;
    y0 = sum(sum(z.*py))/sz;
    disp('gaussfit2DMLE out of bounds: returning (x0, y0) = centroid!')
end

end

    function negL = objfun(params,px,py,z)
    %  Log-likelihood for Poisson-distributed intensity values, with a
    %  Gaussian spatial distribution
        dr = [px(:) - params(2),py(:)-params(3)];
        gaussprob = params(5)*exp(-sum(dr.*dr,2)/2/params(4)/params(4)) + params(1);
        Lk = (z(:)).*log(gaussprob) - gaussprob;  % log-likelihood; Poisson-distr.
        negL = -sum(Lk(:));
    end
    