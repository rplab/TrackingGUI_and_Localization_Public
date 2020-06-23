% cullTracks_function.m
% 
% Function to "edit" tracks based on various criteria.
% Can also simply evaluate, return criteria, without culling.
%
% Replaces trackedit.m, to more easily allow a variety of criteria,
%    addition of criteria, and calling from cullTracksGUI.m
%
% Inputs
%   [Required]
%   objs : linked object matrix, for example from nnlink_rp.m  (See that file, or
%          im2obj_rp for object structure). Typically this is objs_link
%   cullString : string that specifies what track property to evaluate. 
%          NOTE: Case sensitive, but I also include lowercase variants
%          in the case expression
%          Options:
%          'StdDev' : standard deviation of position / sqrt(track length)
%                    (esp. useful for stuck objects). 
%          'TrackLength' : track length (frames)
%          'Straightness' : straightness index (0, 1] -- the ratio of the
%                     end-to-end distance to the contour length
%          'Angle' : the average angle of a track, [0, 360 degrees], 
%                     from best-fit line. 
%          'StepSpeed' : Average frame-to-frame displacement of a track.
%   [Optional -- input as 'Parameter' , Value pairs]
%   perform_culling : if true (default) cull tracks; if false, will still
%           evaluate track properties, and return this
%   dispopt : 2-element array. 
%           1) Show track information (text) if true (optional; default false)
%           2) Show histogram of property values, to help user with criteria
%              Forced to true if params is empty or not input
%   cullParams : culling parameters for the selected property; Leave empty to
%           allow user to input, based on histograms.
%           2-element array; typically [min max] values of the property.
%           For angles, if max < min, assume we're considering a range of
%           angles that crosses 0 (360) degrees, so keep >max OR <min
%
% Output
%   objs_out : object matrix of tracks that are kept
%           Doesn't alter track IDs, so output IDs may not be sequential
%   track_properties : values of the evaluated track property for the
%           culled tracks (if perform_culling is TRUE), or all tracks
%           (otherwise). In the same order as track_IDs.
%   track_IDs : Track IDs (row 6 of object matrix) of the tracks that are
%           retained. This is just unique(objs_out(6,:)), but I'll output
%           it separately anyway.
%   h_hist : handle to histogram figure window.
%   outputString : text string summarizing number of input, output tracks,
%           and parameters; possibly useful for GUI. Display if dispopt(1)
%           is true.
%
% To do
%   - Check screen size, to make sure figure isn't too big
%   - (cosmetic): save binEdges from the subplots (allow for 1 or 2) so that
%     reploting the histogram has the same bins
%
% Raghuveer Parthasarathy
% March 23, 2020 -- based on trackedit.m, and trackveldist.m
% May 8, 2020 -- use Input Parser class to validate inputs
% last modified May 8, 2020

function [objs_out, track_properties, track_IDs, h_hist, outputString] = ...
    cullTracks_function(objs, cullString, varargin)


h_hist = []; % remains empty, if no figure is made

%% Defaults

default_Perform_culling = true;
default_dispopt = [false false];
default_cullParams = [];

% Input Parser
p = inputParser;
validObjsSize = @(x) size(x,1)>=6 && size(x,2)>0; % at least 6 rows and 1 col
validCullStrings = {'StdDev', 'TrackLength', 'Straightness', 'Angle', 'StepSpeed'};
addRequired(p,'objs',validObjsSize);
addRequired(p,'cullString', @(x) any(validatestring(x,validCullStrings)));
addOptional(p,'perform_culling',default_Perform_culling, @(x) islogical(x) || x==1 || x==0);
addOptional(p,'dispopt',default_dispopt,@(x) (islogical(x) || x==1 || x==0) && length(x)<=2);
addOptional(p,'cullParams',default_cullParams,@(x) length(x)<=2);
parse(p, objs, cullString, varargin{:});

% modifications to input parameters
cullString = lower(p.Results.cullString); % make lowercase
dispopt = p.Results.dispopt; % 
cullParams = p.Results.cullParams; % 
if length(dispopt)==1
    dispopt(2) = false;
end
if ~exist('cullParams', 'var') || isempty(cullParams)
    dispopt(2) = true;
end


%% Evaluate criterion

% Different options for different criteria

switch cullString
    case {'stddev'}
        % Calculate standard deviation per frame (px)
        track_property_function = @get_track_StdDev; % handle to function
        xlabelString = 'Std. Dev. / sqrt(Nframes), px'; % for plots
        cullParamDefaults = [0 Inf]; % default min, max for culling range
        
    case {'tracklength'}
        % Calculate length of each track (number of frames)
        track_property_function = @get_trackLength; % handle to function
        xlabelString = 'Track Length, frames';
        cullParamDefaults = [0 Inf]; % default min, max for culling range

    case {'straightness'}
        % Calculate straigtness index of each track 
        track_property_function = @get_Straightness; % handle to function
        xlabelString = 'Straightness Index';
        cullParamDefaults = [0 Inf]; % default min, max for culling range

    case {'angle'}
        % Calculate angle [0, 2*pi) of each track 
        track_property_function = @get_Angle; % handle to function
        xlabelString = 'Angle, degrees';
        cullParamDefaults = [0 360]; % default min, max for culling range

    case {'stepspeed'}
        % Calculate the mean displacement per frame (i.e. speed) of each
        % track
        track_property_function = @get_StepSpeed; % handle to function
        xlabelString = 'Step Speed, px/frame';
        cullParamDefaults = [0 Inf]; % default min, max for culling range
end

%% Evaluate track properties -- call various functions
[track_IDs, ~, iUtrk] = unique(objs(6,:)); % Unique Track IDs, and indexing
   % thisTrackNo = iUtrk(j); % the "unique" track index in utrk
   % corresponding to the object in the j'th column of objs
track_properties = track_property_function(objs, track_IDs, iUtrk);

%% Get culling parameters, if not input
if ~exist('cullParams', 'var') || isempty(cullParams)
    [cullParams, h_hist] = get_culling_params(track_properties, cullString, xlabelString, cullParamDefaults);
else
    if dispopt(2)
        % parameters already input, but show histogram anyway
        h_hist = make_histogram(track_properties, cullString, xlabelString);
    end
end

%% Cull, based on parameters

if p.Results.perform_culling
    % Cull, based on params
    objs_out = cullTracks(objs, track_properties, cullParams);
    [track_IDs, ~, iUtrk] = unique(objs_out(6,:)); % Unique Track IDs, and indexing
    track_properties = track_property_function(objs_out, track_IDs, iUtrk);
    
    % Display properties of retained objects, if desired
    if dispopt(2)
        % Recalculate and Replot
        make_histogram(track_properties, cullString, xlabelString, h_hist);
    end
else
    objs_out = objs;
end
track_IDs = unique(objs_out(6,:));


%%

Ntracks_original = length(unique(objs(6,:)));
% Output string. Use {' '} to preserve spaces in strcat. Then make string
outputString = strcat('culled:', {' '}, cullString, ', Params', {' '}, ...
    sprintf('[%.2f %.2f]; ', cullParams(1), cullParams(2)), {' '}, ...
    sprintf('Input %d tracks; %d rejected; %d kept', ...
       Ntracks_original, Ntracks_original- length(track_IDs), length(track_IDs)));
outputString = char(outputString);

if dispopt(1)
    % Number of starting, ending tracks
    disp(outputString);
    fprintf('Input object array: %d unique tracks.\n', Ntracks_original);
    fprintf('Output object array: %d unique tracks.\n', length(track_IDs));
end


end

%% Functions for evaluating properties

function trk_StdDevPerSqrtFrame = get_track_StdDev(objs, utrk, iUtrk)
% Calculate standard deviation / Number of frames for each track
% Do this by incrementing sum(x) and sum(x^2) arrays, using this to
% calculate variance. This is *far* faster (~100x) than extracting the
% columns corresponding to each track and calculating the variance of each
% set of positions. See "deleted_from_trackedit.m" file -- Mar 11, 2020
    Nframes = zeros(1,length(utrk));
    sum_x = zeros(1,length(utrk));
    sum_x2 = zeros(1,length(utrk));
    sum_y = zeros(1,length(utrk));
    sum_y2 = zeros(1,length(utrk));
    progtitle = sprintf('cullTracks_function: Calculating variance...  ');
    progbar = waitbar(0, progtitle);  % will display progress
    for j=1:size(objs,2)
        % loop through each column
        % Increment sum(x), sum(x2), Nframes for the track this column is
        % part of.
        thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
        % Increment positions; don't keep track of NaNs
        sum_x(thisTrackNo) = sum_x(thisTrackNo) + objs(1,j);
        sum_x2(thisTrackNo) = sum_x2(thisTrackNo) + (objs(1,j))^2;
        sum_y(thisTrackNo) = sum_y(thisTrackNo) + objs(2,j);
        sum_y2(thisTrackNo) = sum_y2(thisTrackNo) + (objs(2,j))^2;
        Nframes(thisTrackNo) = Nframes(thisTrackNo) + 1;
        if (mod(j,round(size(objs,2)/25))==0)
           waitbar(j/size(objs,2), progbar, strcat(progtitle, sprintf('Object %d of %d', j, size(objs,2))));
        end
    end
    close(progbar)
    % Calculate variance (Normalize by N rather than N-1, just to be
    % concise.)
    vx = sum_x2./Nframes - (sum_x./Nframes).^2;
    vy = sum_y2./Nframes - (sum_y./Nframes).^2;
    trk_StdDevPerSqrtFrame = sqrt((vx+vy)./Nframes); 
end

function trk_Length = get_trackLength(objs, utrk, iUtrk)
    % Determine the length of each track. 
    % Avoid "Nframes = sum(objs(6,:)==j);" -- very slow for large arrays
    trk_Length = zeros(1,length(utrk));
    for j=1:size(objs,2)
        thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
        trk_Length(thisTrackNo) = trk_Length(thisTrackNo) + 1;
    end
end
    
function trk_Straightness = get_Straightness(objs, utrk, iUtrk)
    % Determine the Straightness Index of each track. (0, 1]
    % Note: will be NaN for a track that's all zeros
    % As with the standard deviation calculation, increment arrays rather
    % than extracting all columns corresponding to a particular track, for
    % speed. We need to sum dx^2 and dy^2, which requires keeping the
    % previous frame's position.
    % Requires that frame number be increasing with column number for 
    % a given track; hard to imagine that this wouldn't be the case, 
    % but check anyway:
    if min(diff(objs(5,:)) & diff(objs(6,:)==0)) < 0
        % second condition checks that we're looking at the same track
        errordlg('get_Straightness in cullTracks_function: Frame number not increasing with columns in objs!')
    end
    % Avoid "Nframes = sum(objs(6,:)==j);" -- very slow for large arrays

    Nframes = zeros(1,length(utrk));
    start_pos = zeros(2, length(utrk)); % x; y
    end_pos = zeros(2, length(utrk)); % x; y
    prev_pos = zeros(2, length(utrk));
    sum_ds = zeros(1,length(utrk));
    progtitle = sprintf('cullTracks_function: Calculating straightness index...  ');
    progbar = waitbar(0, progtitle);  % will display progress
    for j=1:size(objs,2)
        % loop through each column
        % Increment sum(dx^2),..., by keeping track of previous position
        thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
        if Nframes(thisTrackNo)==0
            % first position
            start_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)];
            % don't calculate dx, dy
        else
            % not the first position
            % Increment contour length; don't keep track of NaNs
            sum_ds(thisTrackNo) = sum_ds(thisTrackNo) + ...
                sqrt((objs(1,j) - prev_pos(1, thisTrackNo)).^2 + ...
                     (objs(2,j) - prev_pos(2, thisTrackNo)).^2);
        end
        end_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)]; % keep overwriting this
        prev_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)];
        Nframes(thisTrackNo) = Nframes(thisTrackNo) + 1;
        if (mod(j,round(size(objs,2)/25))==0)
           waitbar(j/size(objs,2), progbar, strcat(progtitle, sprintf('Object %d of %d', j, size(objs,2))));
        end
    end
    close(progbar)
    % Calculate Straightness Index
    % Could compare contour length to length of best-fit line, but this 
    %   allows (rarely) a straightness index > 1.0. Instead use
    %   start-to-end distance -- noisy, but in (0, 1].
    % [~, ~, B, ~] = fitline(x, y, [], false); % fit trajectory to a line
    % length_line = sqrt(1 + B^2)*(max(x)-min(x)); % length of best-fit line
    length_line = sqrt(sum((end_pos - start_pos).^2));
    trk_Straightness = length_line./sum_ds; 
end

function trk_Angle = get_Angle(objs, utrk, iUtrk)
    % calculate the average angle of a track, 0-360 *degrees*, from best-fit 
    % line to the trajectory
    % As with the standard deviation calculation, increment arrays rather
    % than extracting all columns corresponding to a particular track, for
    % speed. We need to sum x, y, etc., for the linear regression
    % Need to keep track of start and end position, to identify quadrant
    %
    % Requires that frame number be increasing with column number for 
    % a given track; hard to imagine that this wouldn't be the case, 
    % but check anyway:
    if min(diff(objs(5,:)) & diff(objs(6,:)==0)) < 0
        % second condition checks that we're looking at the same track
        errordlg('get_Angle in cullTracks_function: Frame number not increasing with columns in objs!')
    end
    % Avoid "Nframes = sum(objs(6,:)==j);" -- very slow for large arrays

    start_pos = zeros(2, length(utrk)); % x; y
    end_pos = zeros(2, length(utrk)); % x; y
    Nframes = zeros(1,length(utrk));
    sum_x = zeros(1, length(utrk));
    sum_x2 = zeros(1, length(utrk));
    sum_y = zeros(1, length(utrk));
    sum_xy = zeros(1, length(utrk));
    progtitle = sprintf('cullTracks_function: Calculating Angle...  ');
    progbar = waitbar(0, progtitle);  % will display progress
    for j=1:size(objs,2)
        % loop through each column, Increment sums
        thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
        if Nframes(thisTrackNo)==0
            % first position
            start_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)];
        end
        end_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)]; % keep overwriting this
        sum_x(thisTrackNo) = sum_x(thisTrackNo) + objs(1,j);
        sum_x2(thisTrackNo) = sum_x2(thisTrackNo) + objs(1,j)^2;
        sum_y(thisTrackNo) = sum_y(thisTrackNo) + objs(2,j);
        sum_xy(thisTrackNo) = sum_xy(thisTrackNo) + objs(1,j)*objs(2,j);
        Nframes(thisTrackNo) = Nframes(thisTrackNo) + 1;
        if (mod(j,round(size(objs,2)/25))==0)
           waitbar(j/size(objs,2), progbar, strcat(progtitle, sprintf('Object %d of %d', j, size(objs,2))));
        end
    end
    close(progbar)
    slope = (Nframes.*sum_xy - sum_x.*sum_y)./ (Nframes.*sum_x2 - sum_x.^2); % linear regression
    radian_Angle = atan2(slope.*sign(end_pos(2,:)-start_pos(2,:)), ...
        sign(end_pos(1,:)-start_pos(1,:)));  % [-pi, pi]
    radian_Angle(radian_Angle<0) = radian_Angle(radian_Angle<0)+2*pi;
    trk_Angle = radian_Angle*180/pi;
end


function trk_StepSpeed = get_StepSpeed(objs, utrk, iUtrk)
    % Determine the average displacement per frame of each track
    % As with the standard deviation calculation, increment arrays rather
    % than extracting all columns corresponding to a particular track, for
    % speed. We need to sum dx^2 and dy^2, which requires keeping the
    % previous frame's position.
    % Requires that frame number be increasing with column number for 
    % a given track; hard to imagine that this wouldn't be the case, 
    % but check anyway:
    if min(diff(objs(5,:)) & diff(objs(6,:)==0)) < 0
        % second condition checks that we're looking at the same track
        errordlg('get_StepSpeed in cullTracks_function: Frame number not increasing with columns in objs!')
    end
    % Avoid "Nframes = sum(objs(6,:)==j);" -- very slow for large arrays

    Nframes = zeros(1,length(utrk));
    prev_pos = zeros(2, length(utrk));
    sum_ds = zeros(1,length(utrk));
    progtitle = sprintf('cullTracks_function: Calculating straightness index...  ');
    progbar = waitbar(0, progtitle);  % will display progress
    for j=1:size(objs,2)
        % loop through each column
        % Increment sum(dx^2),..., by keeping track of previous position
        thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
        if Nframes(thisTrackNo)>0
            % not the first position
            % Increment contour length; don't keep track of NaNs
            sum_ds(thisTrackNo) = sum_ds(thisTrackNo) + ...
                sqrt((objs(1,j) - prev_pos(1, thisTrackNo)).^2 + ...
                     (objs(2,j) - prev_pos(2, thisTrackNo)).^2);
        end
        prev_pos(:, thisTrackNo) = [objs(1,j); objs(2,j)];
        Nframes(thisTrackNo) = Nframes(thisTrackNo) + 1;
        if (mod(j,round(size(objs,2)/25))==0)
           waitbar(j/size(objs,2), progbar, strcat(progtitle, sprintf('Object %d of %d', j, size(objs,2))));
        end
    end
    close(progbar)
    trk_StepSpeed = sum_ds./Nframes; 
end

function [cullParams, h_hist] = get_culling_params(trk_property, cullString, xlabelString, cullParamDefaults)
% Get parameters for culling, from user input
% Inputs ...
%   ... paramDefaultRange : default range for min, max of parameter
    h_hist = make_histogram(trk_property, cullString, xlabelString);
    % Dialog box for parameters
    prompt = {'Min. value:', 'Max. value (Use max<min for Angles that cross 360 deg.):'};
    dlg_title = strcat('Parameters : ', cullString); num_lines= 1;
    def     = string(cullParamDefaults) ;  % default values   % {'0', 'Inf'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    cullParams = [str2double(char(answer(1))) str2double(char(answer(2)))];
end

function h_hist = make_histogram(trk_property, cullString, xlabelString, h_hist)
% Make a histogram of a track property
% Can input figure handle, if made previously (to update)
% If all values are positive, make both linear and log-spaced hisograms
    if ~exist('xlabelString', 'var') || isempty(xlabelString)
        xlabelString = cullString;
    end
    if ~exist('h_hist', 'var') || isempty(h_hist)
        h_hist = figure('name', strcat('cullTracks_function : ', cullString), ...
        'position', [50 50 900 550]); 
        histogram_Facecolor = [0.3 0.7 0.9];
    else
        figure(h_hist);
        histogram_Facecolor = [0.9 0.6 0.2];
    end
    if (min(trk_property) >= 0) && (max(trk_property) >= 0)
        % All the same sign, positive
        % Make log and linear histogram
        % Note some values of std. dev. will be zero, but I still want to
        % make a log-scale histogram
        subplot(1, 2, 1); % for linear histogram
        histogram(trk_property, 'FaceColor', histogram_Facecolor)
        xlabel(xlabelString);
        ylabel('Number of Tracks');
        hold on
        subplot(1, 2, 2); % for log histogram
        % Adjust bins, in case there are zeros
        minx = floor(log10(min(trk_property(sign(trk_property)==1))));
        maxx = ceil(log10(max(trk_property)));
        Nbins = max(20, (maxx-minx)/3); % number of bins; at least 20
        x=logspace(minx, maxx, Nbins); % create bin edges with logarithmic scale
        histogram(trk_property, x, 'FaceColor', histogram_Facecolor);
        set(gca,'xscale','log'); % scale the x-axis
        xlabel(xlabelString);
        ylabel('Number of Tracks');
        hold on
    else
        histogram(trk_property, 'FaceColor', histogram_Facecolor)
        xlabel(xlabelString);
        ylabel('Number of Tracks')
        hold on
    end
    hold on
end


function objs_out = cullTracks(objs, trackProperties, params)
% Remove tracks with parameter values outside the params range (min, max)
% If max < min, assume that these are Angles, with a desired range to keep
% that crosses 360 degrees.
    utrk = unique(objs(6,:));  % all the unique track ids
    objs_out = zeros(size(objs));  % the largest it could possibly be
    
    nc = 1;  % number of columns, for re-sizing the array
    progtitle = sprintf('Culling tracks...  ');
    progbar = waitbar(0, progtitle);  % will display progress
    for j=1:length(utrk)
        if params(2) >= params(1)
            % Usual [min max]
            keep_track_condition = (trackProperties(j) >= params(1)) && (trackProperties(j) <= params(2));
        else
            % max < min, so assume need to wrap angles around 360
            keep_track_condition = (trackProperties(j) <= params(2)) || (trackProperties(j) >= params(1));
        end
        if keep_track_condition
            % keep this track
            trtmp = objs(:, objs(6,:)==utrk(j));  % objects that are part of track j
            objs_out(:,nc:(nc+size(trtmp,2))-1) = trtmp; % keep these
            %k = k+1;
            nc = nc + size(trtmp,2);
        end
        if (mod(j,round(length(utrk)/25))==0)
            waitbar(j/length(utrk), progbar, strcat(progtitle, sprintf('track %d of %d', j, length(utrk))));
        end
    end        
    close(progbar)
    nc = nc-1;
    objs_out = objs_out(:,1:nc);  % "re-sizing" the array
end

