% fitline.m
%
% function to fit a line y = A + Bx to data sets [x], [y] 
% Returns A, B, and uncertainties.  Also returns the coefficient of determination 
%   (R2), Pearson product-moment correlation coefficient (r), and handle to (optional) figure
% 
% If input, uncertainties in y are used in the determination of A, B (Aug.
% 6, 2009).  If these are not input, uncertainties in A and B are 
% calculated using uncertainties in y as estimated from the deviation
% from the line 
% See Bevington (1969); uncertainties in A, B from Eq. 6-25
%
% Inputs:
%     x : x array
%     y : y array
%     sigy : array of uncertainties in y; optional
%     plotopt : IF plotopt==1, plot the points together with 
%          the best-fit line; default false.
%          Also plots the 1-sigma confidence interval about the best fit
%          line, for 100 points between min and max of x
%          (doesn't output anything about this -- see formula below;
%          can easily do this externally)
%     varargin : variable input arguments for the plot, 
%          e.g. "'marker','o', 'color', 'k', 'markersize', 12"
%          Defaults: 'ko', 'markerfacecolor', [0.2 0.6 1.0], 'markersize', 10
% 
% Outputs
%     A, sigA : Intercept, and uncertainty
%     B, sigB : Slope, and uncertainty
%     R2 : coefficient of determination
%     r  : Pearson product-moment correlation coefficient
%     h  : handle to figure (empty if plotopt==false)
%     If there are <5 output arguments, don't bother calculating r or R2
%          (for speed).  Annoying: if output h is used, function will
%          always calculate r or R2.  Could fix with argument check
%
% Raghuveer Parthasarathy
% modified Nov. 5, 2008 to calc. correlation coefficient (R2)
%          Nov. 23, 2010 (faster; calc. R2 only if desired)
%          Sept 16, 2014 (require option 3 to be sigy)
% modified Feb. 21, 2015 (calc. Pearson's r)
% modified Mar. 27, 2015 (1 std. standard deviation about the best-fit line)
% last modified Dec. 13, 2016 (varargin option)

function [A, sigA, B, sigB, R2, r, h] = fitline(x, y, sigy, plotopt, varargin)


if (length(x(:)) ~= length(y(:)))
    disp('Error (fitline.m)!  x, y are not the same length!')
    input('Recommend Control-C to abort. [Enter]');
end

if ~exist('sigy', 'var') 
    sigy = []; 
end

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end

if ~exist('varargin', 'var') || isempty(varargin)
    varargin = {'ko', 'markerfacecolor', [0.2 0.6 1.0], 'markersize', 10};
end

%%
sigyinput = true;  % is the function called with uncertainty values? Default: true
if (nargin < 4)
    plotopt = false;
    if (nargin < 3)
       sigy = ones(size(x));  % irrelevant, uniform uncertainty
       sigyinput = false;
    else
       % There are three input arguments
       if (length(sigy(:))==1)
          % three arguments, and the third is just one number -- so this argument
          % is probably "plotopt", and the function is being called by something
          % written before the consideration of sigy
          plotopt = sigy;
          sigy = ones(size(x));
          sigyinput = false;
       end
    end
end
if isempty(sigy)
    % in case it's empty
    sigy = ones(size(x));
    sigyinput = false;
end
%%

% To ensure arrays are same shape
x = x(:);
y = y(:);
sigy = sigy(:);
N = length(x);

if (N<2)
    % Zero or 1 elements in array
    A = NaN; sigA = NaN; B = NaN; sigB = NaN; R2 = NaN; h = [];
else
    % Least squares linear fit: y = A + Bx
    sx = sum(x./sigy./sigy);
    sxx = sum(x.*x./sigy./sigy);
    sy = sum(y./sigy./sigy);
    sxy = sum(x.*y./sigy./sigy);
    ssig = sum(1./sigy./sigy);
    D = ssig*sxx - (sx*sx);
    A = (sxx*sy - sx*sxy)/D;
    B = (ssig*sxy - sx*sy)/D;
    if nargout>=5
        % user wants R2  and r to be calculated; the slowest step in the function,
        % so do only if requested.
        % Annoyance: will also do if h is output, since nargout==7 in this
        % case.
        meanx = sum(x)/N;  % faster than mean(x);
        meany = sum(y)/N;  % faster than mean(y);
        dxdy = sum((x-meanx).*(y-meany));
        ddx = sum((x-meanx).^2);
        ddy = sum((y-meany).^2);
        r = dxdy / sqrt(ddx)/sqrt(ddy);
        R2num1 = sum(x.*y) - N*meanx*meany;
        R2 = R2num1*R2num1 / ...
            ((sum(x.*x) - N*meanx*meanx)*(sum(y.*y) - N*meany*meany));
    else
        r = NaN;
        R2 = NaN;
    end
    % Uncertainties
    if (N > 2)
        if (sigyinput == false)
            % use deviation from the line as the uncertainty
            devfromline = y - A - B*x;
            sigyfitsq = sum(devfromline.*devfromline)/(N-2);
            sigA = sqrt(sigyfitsq*sxx/D);
            sigB = sqrt(sigyfitsq*N/D);
        else
            sigA = sqrt(sxx/D);
            sigB = sqrt(ssig/D);
        end
    else
        sigA = NaN;
        sigB = NaN;
    end
    
    if plotopt
        h = figure; plot(x, A + B*x, '-', 'Color', [0.5 0.5 0.5]);
        hold on;
        plot(x, y, varargin{:});
    else
        h = [];
    end

end

% Calculate standard deviation about the best-fit line (see 
% http://www.cs.wayne.edu/~hzhang/courses/7290/Lectures/11%20-%20Simple%20Linear%20Regression%20Models.pdf
% saved in \Statistics
if N>2 && plotopt
    % redundant with above, if sigyinput==false 
    devfromline = y - A - B*x;
    sigyfitsq = sum(devfromline.*devfromline)/(N-2);
    xp = linspace(min(x), max(x), 100);
    yp = A + B*xp;
    mx = mean(x);
    syp = sqrt(sigyfitsq*(1/N + (xp-mx).^2 / (sxx - N*mx*mx)));
    plot(xp, yp+syp, ':', 'Color', 0.7*[1 1 1]);
    plot(xp, yp-syp, ':', 'Color', 0.7*[1 1 1]);
end