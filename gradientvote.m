% gradientvote.m
% 
% function to calculate the gradient at each point of an image, and then
% assign to each pixel along the gradient line within some distance a "vote"
% proportional to the gradient magnitude.  (The gradient line is the line
% of slope = intensity gradient, intersecting the point.)
% Application: transform an image of rings, e.g, into an image of bright
% points, without the parameter-dependence of a Hough transformation
% See notes July 10-13, 2012
%
% Inputs
%   I        : image
%   objsize  : object size, px.  (Shouldn't be smaller than the object;
%              could be larger).  Gradient votes are counted along the 
%              gradient line +/- objsize pixels from each point if 
%              graddir=0; only in one direction if graddir = +/- 1.
%   graddir  : option for the relevant direction of the intensity gradient
%              0 : [default] gradient vote in both directions
%              +1 or -1 : count votes only for pixels in the direction of
%              positive or negative intensity gradient, respectively.
%   grthresh : (optional) gradient magnitude threshold [0-1).  If >0, allow
%              only points in the >=grthresh fraction of gradient
%              magnitudes to "vote".  (Default 0)
%
% Outputs
%   I_grvote : the total gradient vote at each pixel.
%
% Raghuveer Parthasarathy
% July 10, 2012
% last modified July 13, 2012

function I_grvote = gradientvote(I, objsize, graddir, grthresh)


%% Calculate the gradient line at each grid midpoint 
% -- similar calculation as in radialcenter.m, but grid midpoints are 
% [1.5 2.5 3.5 ... Nx-0.5] instead of being centered at 0

if ~exist('graddir', 'var') || isempty(graddir)
    graddir = 0;
end
if ~exist('grthresh', 'var') || isempty(grthresh)
    grthresh = 0;
end

% Make sure I is double to make m double, etc.
I = double(I);

% Number of grid points
[Ny Nx] = size(I);

% grid midpoint coordinates ;
xm_onerow = 1.5:(Nx-0.5);
ym_onecol = (1.5:(Ny-0.5))';  % Note that y increases "downward"

% Calculate derivatives along 45-degree shifted coordinates (u and v)
% Note that y increases "downward" (increasing row number) -- we'll deal
% with this when calculating "m" below.
dIdu = I(1:Ny-1,2:Nx)-I(2:Ny,1:Nx-1);
dIdv = I(1:Ny-1,1:Nx-1)-I(2:Ny,2:Nx);

% Don't smooth
dImag2 = dIdu.*dIdu + dIdv.*dIdv; % gradient magnitude, squared

% Slope of the gradient .  Note that we need a 45 degree rotation of 
% the u,v components to express the slope in the x-y coordinate system.
% The negative sign "flips" the array to account for y increasing
% "downward"
m = -(dIdv + dIdu) ./ (dIdu-dIdv); 

infslope = 9e9;  % replace infinite slope values with this extremely large number
m(isinf(m)) = infslope;

% if the intensity is completely flat, m will be NaN (0/0)
% give these points zero votes, by setting dImag2 = 0 (and set m= 0 to avoid 0*NaN=NaN)
dImag2(isnan(m))=0;
m(isnan(m))=0;

I_grvote = zeros(size(I,1)+2*objsize, size(I,2)+2*objsize);

% create a grid of rj values, in a small array of size 2*(objsize+1) * 2*(objsize+1)
xj_onerow = -(objsize+0.5):(objsize+0.5);
xj = xj_onerow(ones(2*(objsize+1), 1), :);
yj_onecol = (-(objsize+0.5):(objsize+0.5))';  % Note that y increases "downward"
yj = yj_onecol(:,ones(2*(objsize+1),1));
rj = sqrt(xj.*xj + yj.*yj);  % distance to center of the small array
toolargerj = rj>objsize;

% gradient intensity thresholding (only consider points in the top 
% fraction of gradient magnitudes; e.g. if grthresh = 0.9, consider the 
% "brightest" 10% of points)
% First determine indices of above-threshold points
if grthresh > 0.0
    % Do thresholding
    [hs, bins] = hist(dImag2(:),1000);
    ch = cumsum(hs);
    ch = ch/max(ch);
    noiseind = find(ch > grthresh); %
    noiseind = noiseind(1); % The index value below which "grthresh" fraction
    % of the values lie.
    % find indices of dImag2's that are above-threshold
    [goodk,goodj] = find(dImag2 > bins(noiseind));  
    % Note that goodj(3),goodk(3) is the coordinate of above threshold pt. 3, for example
else
    % No thresholding
    [goodj, goodk] = meshgrid(1:(Nx-1),1:(Ny-1));
end
        
for ptindex = 1:length(goodj(:))
    j=goodj(ptindex);
    k = goodk(ptindex);
    % loop through each grid midpoint
    % calculate the array of d^2 values for each pixel in the
    % neighborhood of each grid midpoint
    d2 = (yj-m(k,j)*xj).^2/(m(k,j)^2+1);
    % vote strength
    thisvotestrength = sqrt(dImag2(k,j))*(1-d2);
    switch graddir
        case 0
            % don't care about direction of gradient
            thisvotestrength((d2>1) | toolargerj) = 0;
        case 1
            % orientation of gradient direction: see July 12, 2012 notes
            costheta = xj.*(dIdu(k,j)-dIdv(k,j)) - yj.*(dIdu(k,j)+dIdv(k,j));
            thisvotestrength((d2>1) | toolargerj | costheta<0) = 0;
        case -1
            % orientation of gradient direction: see July 12, 2012 notes
            costheta = xj.*(dIdu(k,j)-dIdv(k,j)) - yj.*(dIdu(k,j)+dIdv(k,j));
            thisvotestrength((d2>1) | toolargerj | costheta>0) = 0;
        otherwise
            disp('** gradientvote.m: bad graddir!');
            pause(2)
            thisvotestrength = [];
    end
    yrows = ym_onecol(k)-0.5:ym_onecol(k)+(2*objsize+0.5);
    xcols = xm_onerow(j)-0.5:xm_onerow(j)+(2*objsize+0.5);
    I_grvote(yrows, xcols) = I_grvote(yrows, xcols) + thisvotestrength;
end

% take the central part of the vote-counting array
I_grvote = I_grvote(objsize+1:Ny+objsize, objsize+1:Nx+objsize);

