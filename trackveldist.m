% trackveldist.m
% 
% function to calculate the distribution of object (track) velocities,
% where velocities are calculated over running averages over some number of
% frames.  
% Returns the mean velocity of each track, and other measures such as
% straightness. 
% Note: Straightness is also calculated in cullTracks_function.m, in a
% faster way.
% Plots a histogram of the mean track velocity (optional).
%
% calls fitline.m
%
% Inputs
%   objs : object matrix (linked, typically objs_link)
%   binsize : number of frame *steps* over which to fit a line to get the velocity
%             if 1, equivalent to considering displacement over adjacent
%             frames (2-frame segments); default 1
%   displayopt : if true, make histogram (default false); also display text
%                output
%   vscale : velocity scale (um/px * frames/sec to convert px/fr to um/s)
%           (default 1) -- used for histogram only
%
% Outputs
%   vtr : mean segment-averaged velocity for each track whose length is at least binsize
%         (px/fr, or um/sec if vscale is used)
%   id  : the Track ID of each analyzed track (row 6 of objs)
%   trackLength  : the length of each analyzed track
%   all_segment_v : the velocity of each segment (px/fr, or um/sec if vscale is used)
%   Ntracks : the number of tracks analyzed
%
% Raghuveer Parthasarathy
% January 13, 2012
% June 27, 2016: Several changes:
%     - warn about tracks with non-consecutive frames
%     - calculate the mean velocity of each segment, not just average of
%       each track
%     - segments of size binsize are non-overlapping
% July 10, 2016: output Ntracks, output vtr in um/sec
% March 11, 2020: Analyze straightness of tracks, angle of tracks
% March 25, 2020: Remove analysis of straightness of tracks, angle of
%     tracks -- move to cullTracks_function.m
% last modified May 11, 2020: Quickly written, inelegant export of
%    velocities, track length to Excel file.

function [vtr, id, trackLength, all_segment_v, Ntracks]...
    = trackveldist(objs, binsize, displayopt, vscale)

if ~exist('binsize', 'var') || isempty(binsize)
    binsize = 1;
end
if ~exist('displayopt', 'var') || isempty(displayopt)
    displayopt = false;
end
if ~exist('vscale', 'var') || isempty(vscale)
    vscale = 1;
    % vscale = 0.1625*90;  % from px/fr to um/s
end

unqtrackIDs = unique(objs(6,:)); % get track numbers
trackLength = zeros(size(unqtrackIDs));
vtr = zeros(size(unqtrackIDs));
id = zeros(size(unqtrackIDs));

disp(' ')
all_segment_v = []; % store all segment velocities
for j=1:length(unqtrackIDs)
    % consider each track
    is_trackj = find(objs(6,:) == unqtrackIDs(j)) ; %faster than ismember(objs(6,:),unqtrackIDs(j));
    trackLength(j) = length(is_trackj);  % length of this track
    objs_trj = objs(:, is_trackj);  % object matrix with just track j
    max_frame_gap = max(diff(sort(objs_trj(5,:))));
    if max_frame_gap>1
        fprintf('Warning: Skipping track %d; non-consecutive frame numbers\n', j);
    end
    % Velocity of segments, and average velocity of each track's segments
    if trackLength(j)>=(binsize+1) && max_frame_gap==1
        [vtr(j), v_eachSegment] = getTrackVelocity(objs_trj, trackLength(j), binsize);
        all_segment_v = [all_segment_v v_eachSegment]; % all segment velocities
    else
        vtr(j) = NaN;
    end
    
    % ID of this track
    id(j) = unqtrackIDs(j);
    
end

% Remove tracks that are too short for binning
trackNaN = find(isnan(vtr));
vtr(trackNaN) = [];
id(trackNaN) = [];
trackLength(trackNaN) = [];
Ntracks = length(vtr);

% Convert to "real" units
vtr = vtr*vscale;
all_segment_v = all_segment_v*vscale;


%% Output to Excel file
% Crude; quickly write since it's useful for Cathy
writeExcel = questdlg('Export velocity, track length to Excel?', ...
    'outputExcel', 'yes', 'no', 'no'); % last item is default
if strcmpi(writeExcel, 'yes')
    % Dialog box for file name
    dlg_title = 'File Name for Excel'; num_lines= 1;
    prompt = {'Output filename -- include .xlsx:'};
    def = {'velocities.xlsx'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    ExcelFilename = char(answer(1));
    xlswrite(ExcelFilename,[vtr' trackLength']); % write to file
end

            
%% Make histogram of track properies
% (for tracks with nonzero velocities)
% Also, display average properties (text)
if displayopt
    dv = 0.5*vscale; % 0.25*vscale;
    vshift = 0.0;  % shift bins, to make plot clearer
    
    % Histogram of each segment's velocity
    binv = (dv/2.0:dv:max(all_segment_v)-dv/2.0) + vshift;
    Nsegment = length(all_segment_v(all_segment_v>0));
    Nv_segment = hist(all_segment_v(all_segment_v>0),binv) / Nsegment;
    figure('name', 'histogram of segment velocities', 'position', [100 100 500 500]);
    color_orange = [0.9 0.6 0.2];
    bar(binv, Nv_segment, 'FaceColor', color_orange)
    if vscale ~= 1
        xlabel('Speed (\mum / s)');
    else
        xlabel('Velocity (px/frame)');
    end
    ylabel('Fraction of segments')
    title('Histogram of segment velocities')

    % Histogram of each track's mean velocity
    binv = (dv/2.0:dv:max(vtr)-dv/2.0) + vshift;
    Ntr = length(vtr(vtr>0));
    Nv = hist(vtr(vtr>0),binv) / Ntr;
    figure('name', 'histogram of track velocities', 'position', [500 100 500 500]);
    color_blue = [0.2 0.6 0.9];
    % hold on; color_orange = [0.9 0.6 0.2];
    bar(binv, Nv, 'FaceColor', color_blue)
    if vscale ~= 1
        xlabel('Speed (\mum / s)');
    else
        xlabel('Velocity (px/frame)');
    end
    ylabel('Fraction of trajectories')
    title('Histogram of track velocities')
    
    disp(' '); disp(' ' )
    fprintf('Number of tracks: %d\n', Ntracks); 
    fprintf('Each segment: Mean velocity %.2f +/- %.2f (std.dev.) px/fr or um/s\n', ...
        mean(all_segment_v), std(all_segment_v));
    fprintf('              Median velocity %.2f px/fr or um/s\n', ...
        median(all_segment_v));

    fprintf('Each track: Mean velocity %.2f +/- %.2f (std.dev.) px/fr or um/s\n', ...
        mean(vtr), std(vtr));
    fprintf('            Median velocity %.2f  px/fr or um/s\n', ...
        median(vtr));
end

end

function [v, v_eachSegment] = getTrackVelocity(objs_trj, trackLength, binsize)
% Calculate the velocity for a track
% output the average v of all the segments, and v of each segment
    Nsegments = floor(trackLength/binsize)-1;
    v_eachSegment = zeros(1, Nsegments);
    for k=1:Nsegments % each bin
        x = objs_trj(1,k:k+binsize);  % x positions in this segment of length binsize+1 frames
        y = objs_trj(2,k:k+binsize);  % y positions in this segment of length binsize+1 frames
        if binsize == 1
            % simple frame-to-frame difference
            vx = diff(x);
            vy = diff(y);
        else
            [~, ~, vx, ~] = fitline(1:(binsize+1), x, [], false);
            [~, ~, vy, ~] = fitline(1:(binsize+1), y, [], false);
        end
        v = sqrt(vx*vx+vy*vy);
        v_eachSegment(k) = v;
    end
    v = mean(v_eachSegment); % Average of all segment velocities
end


