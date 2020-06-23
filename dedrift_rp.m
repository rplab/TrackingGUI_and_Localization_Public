function [objs_dedrift, drift] = dedrift_rp(objs, calcMethod, applyMethod, dispopt)
%
% Drift correction for tracked objects.

% Method: determine positions of objects that are in *both*
% of a pair of adjacent frames. Calculate the shift in their positions
% for all image frames, and correct for it (see below).
% Note that this should be robust to gain, loss of particles, and doesn't
% assume any particular form for the drift motion.

% Multiple options for calculating the drift correction (see below)
%
% Inputs
% objs : Object matrix from im2obj_rp or nnlink_rp, with the following form:
%    [x;
%     y;
%     mass;
%     particleid;
%     frame;
%     trackid]
% calcMethod : method for calculating the drift correction
%     Method 1 : Offset by the mean shift of the particles  
%     Method 2 : (default) Offset by the median shift of the particles  
%     Method 3 : Position dependent shift: linear in x, y  
% applyMethod : method for applying the drift correction
%     Method 1 (default) : apply frame-by-frame offset  
%     Method 2 : Average over all frames, then apply  
% dispopt : if true (default), display shift for each frame
%
% Outputs
% objs_dedrift : de-drifted object matrix
% drift : frame-to-frame shift, px, applied (with averaging if applyMethod==2)
%      to each object in each frame 
%      Nframes rows
%      Columns 1 and 2: position independent x and y shifts
%      If calcMethod == 3, 
%         columns 3, 4 are alpha_x (coefficients of x-shift, for x and y terms) 
%         columns 5, 6 are alpha_y (coefficients of y-shift, for x and y terms) 
%
% Raghuveer Parthasarathy
% Major conceptual rewrite, May 24, 2017. See notes.
%   (replaces very simple 2007 version, which just subtracted the median x
%   and y positions from each frame).
% Major changes January, 2019. Allow postition-dependent drift (i.e.
%   shear), and alter outputs.
% Last modified March 22, 2020 (change default to median shift)

%% Defaults
if ~exist('calcMethod', 'var') || isempty(calcMethod)
    calcMethod = 2;
end
if ~exist('applyMethod', 'var') || isempty(applyMethod)
    applyMethod = 1;
end
if ~exist('dispopt', 'var') || isempty(dispopt)
    dispopt = true;
end

%% Calculate drift

% Approximate midpoints, for clarity of referring to the shear terms
% (Could get from the images, but these aren't input
midx = 0.5*max(objs(1,:));
midy = 0.5*max(objs(2,:));

unqframes = sort(unique(objs(5,:))); % all frames
% Check that there are no gaps in frames
if max(diff(unqframes))>1
    warndlg('dedrift_rp.m: WARNING -- gap in frames numbers. Will still run.')
end
delta = zeros(length(unqframes),2); % Constant (no shear) drift at each frame, x and y
alpha_x = zeros(length(unqframes),2); % drift, shear component of x, x and y dependence
alpha_y = zeros(length(unqframes),2); % drift, shear component of y, x and y dependence
progtitle = 'Calculating shift';
progbar = waitbar(0, progtitle);
prevframe_objs = objs(:,objs(5,:)==unqframes(1)); % object matrix of objects in the first frame
for j=2:length(unqframes)
    currframe_objs = objs(:,objs(5,:)==unqframes(j)); % object matrix of objects in the current frame
    currIDs = currframe_objs(6,:); % Track IDs in the current frame
    prevIDs = prevframe_objs(6,:); % Track IDs in the previous frame
    sharedprev = intersect(prevIDs, currIDs); % Track IDs that are in both frames
    % Extract x and y positions in the previous frame of objects that are in the
    % previous and current frames. Caution: ordering of Track IDs may
    % change, so we need to extract each pair (current and former)
    % separately
    x1  = prevframe_objs(1,ismember(prevframe_objs(6,:), sharedprev));
    y1  = prevframe_objs(2,ismember(prevframe_objs(6,:), sharedprev));
    ID1 = prevframe_objs(6,ismember(prevframe_objs(6,:), sharedprev));
    % Extract x and y positions in the current frame of objects that are in the
    % previous and current frames. 
    x2 = zeros(size(x1));
    y2 = zeros(size(y1));
    for k=1:length(x1)
        x2(k) = currframe_objs(1,currframe_objs(6,:)==ID1(k));
        y2(k) = currframe_objs(2,currframe_objs(6,:)==ID1(k));
    end
    
    dx = x2-x1; % shift in x, for these particles
    dy = y2-y1; % shift in y, for these particles

    switch calcMethod
        case 1
            delta(j,:) = [mean(dx), mean(dy)];
        case 2
            % Median, rather than mean
            delta(j,:) = [median(dx), median(dy)];
            % Very similar
        case 3
            % Allow Shear
            % Use initial frame position as the position at which to 
            % calculate shear, for simplicity of correction, below.
            xpos = x1; ypos = y1;
            % Alternatively, could use midpoint.
            % xpos = 0.5*(x2+x1); ypos = 0.5*(y2+y1);
            A = [ones(length(xpos),1) (xpos'-midx) (ypos'-midy)];
            tempSolution = A\dx';
            delta(j,1) = tempSolution(1);
            alpha_x(j,1) = tempSolution(2);
            alpha_x(j,2) = tempSolution(3);
            tempSolution = A\dy';
            delta(j,2) = tempSolution(1);
            alpha_y(j,1) = tempSolution(2);
            alpha_y(j,2) = tempSolution(3);
        otherwise
            errordlg('Error in dedrift_rp.: bad methodopt')
    end
    
    prevframe_objs = currframe_objs;
    if mod(j,20)==0
        waitbar(j/length(unqframes), progbar, progtitle)
    end    
end
close(progbar)

%% Outputs

drift = delta;
if calcMethod==3
    drift = [delta alpha_x alpha_y];
end

%% Apply drift correction

% correct x and y values in the object matrix
objs_dedrift = objs;
% For calculation methods 1 or 2, this is easy -- simple cumulative shift.
% objs_dedrift(5,j) is the frame number of column j; it's probably the same
% as unqframes(objs_dedrift(5,j)), if frame numbers start at 1 and are
% consecutive, but just to be safe let's look for the match
switch applyMethod
    % frame-to-frame shift, or average over frames
    case 1
        TotalShift = cumsum(delta, 1); % cumulative shift at each examined frame
        for j=1:size(objs_dedrift,2)
            calcframe_index = find(unqframes == objs_dedrift(5,j));
            objs_dedrift(1,j) = objs_dedrift(1,j) - TotalShift(calcframe_index,1);
            objs_dedrift(2,j) = objs_dedrift(2,j) - TotalShift(calcframe_index,2);
        end
    case 2
        % Use the average frame-to-frame shift
        meandelta = mean(delta, 1);
        for j=1:size(objs_dedrift,2)
            objs_dedrift(1,j) = objs_dedrift(1,j) - (objs_dedrift(5,j)-min(unqframes))*meandelta(1);
            objs_dedrift(2,j) = objs_dedrift(2,j) - (objs_dedrift(5,j)-min(unqframes))*meandelta(2);
        end
    otherwise
        disp('dedrift_rp.m: Bad applyMethod')
end

if calcMethod==3
    progtitle = 'Applying shift; shear-corrected';
    % For shear, need to consider each trajectory and frame individually, since
    % may move through regions of different motion
    % First, overall (position independent) shift, performed above
    % Then, for each trajectory...
    if applyMethod==2
        % Average over frames
        % Lazily, I'll average and then repeat to use the same indexing as
        % below.
        alpha_x_toUse = repmat(mean(alpha_x,1), length(unqframes), 1);
        alpha_y_toUse = repmat(mean(alpha_y,1), length(unqframes), 1);
    else
        % use frame-to-frame shifts
        alpha_x_toUse = alpha_x;
        alpha_y_toUse = alpha_y;
    end
    unqTracks = unique(objs_dedrift(6,:));
    progbar = waitbar(0, progtitle);
    for j=1:length(unqTracks)
        thisObjs = objs_dedrift(:,objs_dedrift(6,:)==unqTracks(j));
        % Note that thisObjs(5,k) is the frame number; should corresponds to
        % the alphas calculated above for consecutive frames, but make sure
        if size(thisObjs,2)>= 2
            xshift = zeros(1, size(thisObjs,2)); % in each frame, the correction for this track
            yshift = zeros(1, size(thisObjs,2));
            for k=2:size(thisObjs,2)
                calcframe_index = find(unqframes == thisObjs(5,k));
                xshift(k) = alpha_x_toUse(calcframe_index,1)*(thisObjs(1,k-1)-midx) ...
                    + alpha_x_toUse(calcframe_index,2)*(thisObjs(2,k-1)-midy);
                yshift(k) = alpha_y_toUse(calcframe_index,1)*(thisObjs(1,k-1)-midx) ...
                    + alpha_y_toUse(calcframe_index,2)*(thisObjs(2,k-1)-midy);
            end
            thisObjs(1,:) = thisObjs(1,k) + cumsum(xshift);
            thisObjs(2,:) = thisObjs(2,k) + cumsum(yshift);
            objs_dedrift(:,objs_dedrift(6,:)==unqTracks(j)) = thisObjs;
        end
        if mod(j,20)==0
            waitbar(j/length(unqTracks), progbar, progtitle)
        end
    end
    close(progbar)
end

%% Display

if dispopt
    figure;
    plot(unqframes, delta(:,1), '-', 'color', [0.1 0.4 0.7])
    hold on
    plot(unqframes, delta(:,2), '-', 'color', [0.1 0.8 0.3])
    xlabel('Frames')
    ylabel('\delta (px)')
    legend('x', 'y')
    title('Position-independent Drift')
    
    if calcMethod==3
        figure
        plot(unqframes, alpha_x(:,1), '-', 'color', [0.9 0.8 0.4]);
        hold on
        plot(unqframes, alpha_x(:,2), '-', 'color', [0.9 0.6 0.2]);
        legend('\alpha_{xx}', '\alpha_{xy}')
        figure
        plot(unqframes, alpha_y(:,1), '-', 'color', [0.3 0.6 0.9]);
        hold on
        plot(unqframes, alpha_y(:,2), '-', 'color', [0.5 0.7 1.0]);
        legend('\alpha_{yx}', '\alpha_{yy}')
    end
end


