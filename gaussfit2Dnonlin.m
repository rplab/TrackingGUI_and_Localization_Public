% gaussfit2Dnonlin.m
%
% function to fit a 2D (symmetric) gaussian using non-linear regression
%
% based on a post by John D'Errico at
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/165035
%
% fit to form: z = A*exp(-[(x-x0)^2 + (y-y0)^2] / (2*sigma^2)) + offset
% Assume offset >= 0
% Assume same width (sigma) in x, y
%
% Input: 
%    z : 2D array (e.g. particle image)
%    tolerance in z for fitting (default 0.01 if empty)
%              Suggest 0.01 for 8-bit, non-normalized images
%    params0 : [optional] starting values of the parameters
%              If not input or empty, default values
%              1 - constant offset (minimal value of z)

%              2 - x-coordinate of the "center of mass" of the 4 brightest pixels
%              3 - y-coordinate of the "center of mass" of the 4 brightest pixels

%              2 - mode_x  (center of image)
%              3 - mode_y  (center of image)
%              4 - sigma   (quarter of image, roughly)
%              5 - amplitude (max. z - min. z)
%    LB     : [optional] lower bounds of search parameters
%              If not input or empty, default values
%              1 - constant offset (0)
%              2 - mode_x  (px. 1)
%              3 - mode_y  (px. 1)
%              4 - sigma   (0)
%              5 - amplitude (0)
%    UB     : [optional] upper bounds of search parameters
%              If not input or empty, default values
%              1 - constant offset (max value of z)
%              2 - mode_x  (end px.)
%              3 - mode_y  (end px.)
%              4 - sigma   (image size)
%              5 - amplitude (Inf.)
%    lsqoptions : [optional] options structure for nonlinear least-squares
%              fitting, from previosly running 
%              "lsqoptions = optimset('lsqnonlin');"
%              Inputting this speeds up the function.
% Outputs
%    A  : Gaussian amplitude
%    x0, y0 : Gaussian center (in px, from px #1 = left/topmost pixel)
%             So a Gaussian centered in the middle of a 2*N+1 x 2*N+1
%             square (e.g. from make2Dgaussian.m with x0=y0=0) will return
%             a center value at x0=y0=N+1.
%    sigma: Std dev. of Gaussian (assumed same for x, y)
%    offset : constant offset
%
% Raghuveer Parthasarathy 
% March 20, 2011
% July 28, 2011; minor change to default initial parameter values
% Last modified July 28, 2011

function [A, x0, y0, sigma, offset] = gaussfit2Dnonlin(z, tolz, params0, LB, UB, lsqoptions)

[ny,nx] = size(z);
[px,py] = meshgrid(1:nx,1:ny);

% defaults for initial parameter values, and lower and upperbounds
if ~exist('tolz', 'var') || isempty(tolz)
    tolz = 1e-2;
end
if ~exist('params0', 'var') || isempty(params0)
    [B, iB] = sort(z(:), 1, 'descend');
    sB4 = sum(B(iB(1:4)));
    params0_2 = sum(B(iB(1:4)).*px(iB(1:4)))/sB4;
    params0_3 = sum(B(iB(1:4)).*py(iB(1:4)))/sB4;
    params0 = [min(z(:)), params0_2, params0_3, min([nx ny])/4, max(z(:))-min(z(:))];
end
if ~exist('LB', 'var') || isempty(LB)
    LB = [0,1,1,0,0];
end
if ~exist('UB', 'var') || isempty(UB)
    UB = [max(z(:)),nx,ny,max([nx ny]),inf];
end
if ~exist('lsqoptions', 'var') || isempty(lsqoptions)
    lsqoptions = optimset('lsqnonlin');
end

% More fitting options
lsqoptions.TolFun = tolz;  %  % MATLAB default is 1e-6
lsqoptions.TolX = 1e-5';  % default is 1e-6
lsqoptions.Display = 'off'; % 'off' or 'final'; 'iter' for display at each iteration
params = lsqnonlin(@(P) objfun(P,px,py,z),params0,LB,UB,lsqoptions);
A = params(5);
x0 = params(2);
y0 = params(3);
sigma = params(4);
offset = params(1);

end

    function resids = objfun(params,px,py,z)
        temp = [px(:) - params(2),py(:)-params(3)];
        pred = params(1) + params(5)*exp(-sum(temp.*temp,2)/2/params(4)/params(4));
        resids = pred - z(:);
    end
