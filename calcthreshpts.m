% calcthreshpts.m
% 
% Function to determine the points that "pass" the threshold criterion in
% an image.  
% Usual application: particle tracking.  Typically first filter the image 
% (e.g. w/ bpass.m) before inputting to this function.  This function
% thresholds (various options) and looks for local maxima; then, we will look
% at neighborhoods around these with which to 
% localize particles by analyzing the original (non-filtered) image.
%
% Defining the "Local" neighborhood for finding local maxima is a bit
% challenging: large neighborhoods reduce the effect of noise, but can
% obscure dim particles near bright ones.  Choose nsize/4 (one quarter the
% neighborhood size) as a radius of the dilation structuring element (disk),
% with a minimum value of 1.  Previously: was nsize/2. -- RP April 8, 2014
%
% Three thresholding options, that determine the "thresh" variable.
%
% Inputs:
%    A : 2D image.
%    threshopt: determined by fo4_rp.m
%      (1) keeps all pixels with intensity above the thresh*100 percentile
%      (2) keeps pixels with intensity > thresh*std.dev. above
%          the median intensity
%      (3) keeps brightest "thresh" number of particles by finding the 
%          "thresh" brightest local maxima.  Removes local maxima within
%          a neighborhood (radius nsize/4) of other (brighter) maxima.
%          Note that the number of particles can be less than thresh, if
%          fewer than thresh local maxima exist.
%    thresh : intensity threshold
%          For option 2, sign should be positive (flipped by fo5_rp.m) --
%          the function makes sure of this.
%    nsize : necessary only for threshopt==3. object size used for dilating,
%        determining neighborhood size.  Roughly, the particle diameter.
%        Referred to as nsize in > June 27, 2012 fo4_rp.m,
%        since filtering object size is separated from neighborhood object
%        size; formerly called objsize.
%    try1pernhood : (optional; default false).  If true, impose
%        only one local max per neighborhood, by removing local maxima within
%          a neighborhood (radius nsize/2) of other (brighter) maxima.  
%        (This is always done for threshold option 3, regardless of
%        this input variable.)
%    dimg : (optional) dilated image A (pre-calculated for speed)
%          
% Outputs
%    [y, x] : positions of local post-threshold maxima 
%             (N x 1 arrays, where N is the number of maxima.  Note x and y
%             are whole numbers).
%
% Raghuveer Parthasarathy
% extracted from fo4_rp.m, May 2, 2012
% June 28, 2012: write "try1pernhood" option
% February 28, 2013: allow input of dilated image
% March 31, 2014: significant changes to threshold option 3
% April 8, 2014: Alter radius for local maxima finding
% April 9, 2014: Implement new approach to deleting close local maxima, if
%     a brighter one is nearby
% July 18, 2014: For deleting close local maxima, use nsize/2 rather
%       than nsize/4 radius
% Sept. 27, 2019: Avoid error if there are no above-threshold points
%       (options 1 and 2); should return empty x, y arrays
%
% Last modified: Sept. 27, 2019

function [y, x] = calcthreshpts(A, threshopt, thresh, nsize, try1pernhood, dimg)

if ~exist('try1pernhood', 'var') || isempty(try1pernhood)
    try1pernhood = false;
end
if ~exist('dimg', 'var') 
    dimg = [];
end

if threshopt==2
    thresh = abs(thresh);  % make sure it's been flipped
end

showplots=false;  % for debugging -- plot things.

if isempty(dimg)
    % dilate image
    local_disk_radius= max(floor(nsize/4),1);
    ste = strel('disk', local_disk_radius,0);  % for dilation
    dimg = imdilate(A, ste);
end
if showplots
    figure(3)
    imshow(dimg,[]); title('3 dilated')
end
% finding regional maxima (imregionalmax(A)) is much faster than A==dilation.
% However, it is inaccurate, as it does not get rid of local maxima that
% are close (within nsize/2) to brighter true particle maxima.
% imgmax = imregionalmax(A);

% Local maxima in the filtered image
BW = (A == dimg);
if showplots
    figure(4)
    imshow(BW); title('4 local maxima')
end

switch threshopt
    case 1
        % intensity thresholding (percentile)
        [hs, bins] = hist(A(:),1000);
        ch = cumsum(hs);
        ch = ch/max(ch);
        noiseind = find(ch > thresh); %
        noiseind = noiseind(1); % The index value below which "thresh" fraction
        % of the pixels lie.  (Originally noisind(2), but this can lead to errors
        % if there is only one bin above threshold.)
        % find local maxima that are also above-threshold
        % note each local maximum
        [y,x] = find(BW & (A > bins(noiseind)));
        if try1pernhood && ~isempty(x)
            % Delete close maxima
            [x, y] = delete_close_maxima(A, x, y, nsize, []);  % move to function
        end
        if showplots
            figure(5)
            imshow(zeros(size(BW))); hold on; 
            plot(x,y,'wo', 'markersize', 2, 'markerfacecolor', 'w'); 
            title('5 local maxima > threshold')
        end
    case 2
        % intensity thresholding (std. dev. above background)
        medint = median(A(:));
        stdint = std(A(:));
        % find local maxima that are also above-threshold
        isbright = (A > (medint + thresh*stdint));
        % look at each local maximum
        [y,x] = find(BW & isbright);
        if try1pernhood && ~isempty(x)
            % Delete close maxima
            [x, y] = delete_close_maxima(A, x, y, nsize, []);  % move to function
        end
    case 3
        % thresh >= 1 : keep brightest "thresh#" of particles
        % Can't just look for the thresh# of brightest spots -- too sensitive
        % to bright noise.  Need to consider "regions"
        % Method (March 2014): find all local maxima; rank by brightness;
        % for each, starting at the brightest, remove other local
        % maxima closer than nsize/4; repeat
        
        [yy,xx] = find(BW);  % locations of *all* local maxima
        
        % Delete close maxima
        [xx, yy, Nmax] = delete_close_maxima(A, xx, yy, nsize, thresh);  % move to function
        thresh = round(thresh);  % make sure it's an integer
        
        numtofind = min([thresh Nmax]);  % the desired (or total) number of local maxima, and hence objects, to find
        x = xx(1:numtofind);
        y = yy(1:numtofind);
                
end

        function [xx, yy, Nmax] = delete_close_maxima(A, xx, yy, nsize, Nobjs)
        % Method (March 2014): find all local maxima; rank by brightness;
        % for each, starting at the brightest, remove other local
        % maxima closer than nsize/2; repeat
        % Inputs
        %  see main function
        %  Nobjs:  "thresh" for option 3: the maximum number of objects to
        %  "find"  Limit the search to 6*Nobjs, if this is not empty.
        %  Leave empty for options 1, 2
            Amax = A(sub2ind(size(A),yy,xx));  % intensities of all the maxima
            Nmax = length(xx);  % number of maxima found
            [~, fmix] = sort(Amax,'descend');
            % Rank positions based on intensity
            xx = xx(fmix);
            yy = yy(fmix);
            % figure; imagesc(A); colormap('gray'); hold on; plot(xx, yy, 'yo');
            % figure('name', '1'); imagesc(dimg); colormap('gray'); hold on; plot(xx, yy, 'yo');
            
            % Calculate all pairwise distances
            % can be very slow, if there are lots of points!
            % If Nobjs is not empty (option 3), use 6*Nojs as the max
            % extent
            if ~isempty(Nobjs)
                xx = xx(1:min(Nmax,6*Nobjs));
                yy = yy(1:min(Nmax,6*Nobjs));
                Nmax = length(xx);
            end

            % Euclidean distance matrix.
            % uses method of "distance.m" by Roland Bunschoten (MATLAB file exchange)
            xxyy = [xx'; yy']; % should be 2 x Nmax array
            aa=sum(xxyy.*xxyy,1);
            oa = ones(size(aa,2),1);
            % square of distances (matrix)
            if Nmax > 1000
                disp('calcthreshpts.m:  Note that distance matrix calculation can be slow for large numbers of particles.')
            end
            dmat2 = abs(aa( oa, :)' + aa( oa, :) - 2*(xxyy'*xxyy));
            % Inelegant use of a for-loop, and a size-changing array,
            % but avoids improper deletions -- see March 31, 2014 notes
            forNan = tril(ones(size(dmat2)));
            dmat2(forNan==1)=NaN;  % array of NaNs in lower diagonal
            if Nmax > 1
                for j=2:Nmax
                    % examine upper triangle of dmat matrix
                    dmat2(:,dmat2(j,:)<(nsize/2.0)^2)=NaN;  % NaNs in all columns corresponding to a close, dim particle
                    dmat2(dmat2(j,:)<(nsize/2.0)^2,:)=NaN;  % NaNs in all rows corresponding to a close, dim particle
                end
            end
            % delete close particles from the
            yy(isnan([1 dmat2(1,2:end)])) = [];  % get rid of these maxima; always keep particle 1
            xx(isnan([1 dmat2(1,2:end)])) = [];
            % A(isnan([1 dmat(1,2:end)])) = [];
            if ~isempty(Nobjs)
                xx = xx(1:min(length(xx),Nobjs));
                yy = yy(1:min(length(xx),Nobjs));
                Nmax = length(xx);
            end
            % plot(xx, yy, 'rd', 'markerfacecolor', [0.4 1.0 0.7])
        end

end
