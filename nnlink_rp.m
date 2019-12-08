function objs = nnlink_rp(objs, step, memory, dispopt, diagnostic)
% links objects into paths, using nearest-neighbor identification
% with rudimentary checking for unique neighbors via mapping forward &
% backward across frames
% 
% Inputs:
%   objs : object array (see fo5_rp.m for structure; typically output by
%          im2obj_rp.m)
%   step : Cull all links that are greater than sqrt(step*(memory+1)) away
%          Default Inf
%   memory : attempt to link particles despite an absence over some number 
%          of frames, preserving their ID numbers.  Memory==0 (default)
%          means connect adjacent frames only, culling gaps of >= 1 frame.
%          memory==1 allows gaps of 1 frame and culls gaps of >= 2 frames
%   dispopt : display progress if true
%   diagnostic : plots and text to diagnose how this is working
%
% Outputs
%   objs : Object matrix, with row 6 filled in as the track ID
% 
% Raghuveer Parthasarathy
% 16 April, 2007
% May 10, 2007: Fixed foward-backward mapping
% March 16, 2011: Avoid re-numbering frames to start with 1.  
%     Fix treatment of "memory" variable: memory==1 now means consider only
%           linkages spanning 1 frame, culling gaps of 1 frame or more.
% May 21, 2012 -- avoid call to distance.m
% July 18, 2014: Don't remove "lost" objects from the object matrix!
% July 10, 2016: Fixing memory, yet again! (Actually fine, but change 
%      default to zero (no memory), and properly comment)
% last modified July 10, 2016


if ~exist('step', 'var') || isempty(step)
    step = Inf;
end
if ~exist('memory', 'var') || isempty(memory)
    memory = 0;
end
if ~exist('dispopt', 'var') || isempty(dispopt)
    dispopt = false;
end
if ~exist('diagnostic', 'var') || isempty(diagnostic)
    diagnostic = false;
end

unqframes = unique(objs(5,:)); % get frame numbers

startframe = min(unqframes);  % the first frame number (need not be 1)
leftovers = [];

% first assign initial track ids to all objects in the first frame
fm1numobjs = size(objs(4, objs(5,:) == startframe), 2);
objs(6, objs(5,:) == startframe) = 1:fm1numobjs;
nextid = fm1numobjs + 1; % the next available track id

if dispopt
    progtitle = sprintf('nnlink_{rp}: Linking objects...  '); 
    progbar = waitbar(0, progtitle);  % will display progress
end

if diagnostic
    figure('name', 'Initial Track IDs')
    plot(objs(1, objs(5,:) == startframe),objs(2, objs(5,:) == startframe), 'bx')
    for j=1:fm1numobjs
        text(objs(1,objs(6,:) == j), objs(2,objs(6,:) == j), num2str(j));
    end
end

for k = 1:length(unqframes)-1
    cur = [find(objs(5,:) == unqframes(k)) leftovers];
    nex = find(objs(5,:) == unqframes(k+1));     
    % cur and nex are indices of columns in the the object matrix (objs)
    % cur contains all indices corresponding to the present frame 
    % (unqframes(k)), and all the "leftovers" from earlier frames 
    % within the "memory" time
    ncur = length(cur);
    nnex = length(nex);
    if or((ncur == 0),(nnex == 0)) % no elements to track, go to next frame
        disp('No objects to link');
        continue;
    end
    if diagnostic
        disp(' ')
        fs = sprintf('Diagnostic: k = %d, frame %d', k, unqframes(k)); disp(fs)
        fs = sprintf('Diagnostic: Number of objects in the current frame ncur = %d', ncur); disp(fs)
        fs = sprintf('Diagnostic: Number of objects in the next frame nnex = %d', nnex); disp(fs)
    end
    % nearest-neighboor matching algorithm
    pcur = objs(1:2, cur); % positions, current frame
    pnex = objs(1:2, nex); % positions, next frame
    % Euclidean distance matrix.  
    % uses method of "distance.m" by Roland Bunschoten (MATLAB file exchange)
    aa=sum(pnex.*pnex,1); bb=sum(pcur.*pcur,1);
    d = sqrt(abs(aa( ones(size(bb,2),1), :)' + bb( ones(size(aa,2),1), :) - 2*pnex'*pcur));
    % d = distance(pnex,pcur);   % Replace with the above line to avoid extra function call
    % d should be ncur by nnex array of distances
    % d(j,k) is the distance between pnex(j) and pcur(k)
    [~, fmap] = min(d,[],1);  % fmap(j) is the index number in pnex of the particle
                         % that is closest to particle #j in pcur 
                         % fmap has length ncur
                         % Don't use min(d), returns the minimum from each
                         % column of d only if d is not a row vector
    [~, bmap] = min(d,[],2);  % bmap(j) is the index number in pcur of the particle
                         % that is closest to particle #j in pnex 
                         % bmap has length nnex
    invertible = bmap(fmap); % fmap indexes into bmap

    fmapind = 1:length(fmap);
    notlost = fmapind(ismember(fmapind, invertible));  
       % "ismember" is true for a particle of pcur that is mapped back onto
       % itself.  Keep only these in "notlost"
    % Andy Demond's version:
    % invertible(invertible ~= fmapind) = []; % just keep particles for which bmap == inverse(fmap)
    
    if diagnostic
        figure('name', 'mapping')
        plot(pcur(1,:), pcur(2,:), 'kx')
        hold on
        plot(pnex(1,:), pnex(2,:), 'bx')
        plot(pcur(1,invertible), pcur(2,invertible), 'ro')
        title(num2str(k))
    end
    
    % Cull all links that are greater than sqrt(step*(memory+1)) away
    invo1x = objs(1, cur(notlost));
    invo1y = objs(2, cur(notlost));
    invo2x = objs(1, nex(fmap(notlost)));
    invo2y = objs(2, nex(fmap(notlost)));
    sqdisp = (invo1x - invo2x).^2 + (invo1y - invo2y).^2;

    notlost = notlost(sqdisp < step*(memory+1)); %+1 to avoid Inf*0 error.
    if max(size(notlost))==0
        % bizarre, but call this a "1-by-0" empty matrix for the "newguys"
        % calculation below to work
        notlost = zeros(1,0);
    end

%     figure; plot(pcur(2,:), 'ko');
%     hold on; plot(pnex(2,:), 'rs');
%     k
%     fmap
%     bmap'
%     invertible'
%     notlost
%     figure;
%     plot(pcur(1,:), pcur(2,:), 'kx');
%     hold on
%     plot(pnex(1,:), pnex(2,:), 'mx');
%     for jj = 1:length(notlost)
%         plot([pcur(1,notlost(jj)) pnex(1,fmap(notlost(jj)))], [pcur(2,notlost(jj)) pnex(2,fmap(notlost(jj)))],'ko-');
%     end
%     pause
% 
    % assign track ids to particles in next frame
    objs(6, nex(fmap(notlost))) = objs(6, cur(notlost));
    
    if diagnostic
        figure('name', 'Assigning Track IDs')
        plot(objs(1, objs(5,:) == unqframes(k+1)),objs(2, objs(5,:) == unqframes(k+1)), 'yo')
        hold on
        for j=unique(objs(6, cur(notlost)))
            text(objs(1,objs(6,:) == j & objs(5,:) == unqframes(k+1)), ...
                objs(2,objs(6,:) == j & objs(5,:) == unqframes(k+1)), num2str(j));
        end
        title(num2str(unqframes(k+1)))
    end
    
    % cull old objects (i.e. gaps of "memory" frames or more)
    cur(notlost) = [];
    curframes = objs(5, cur);
    cullind = find(curframes < (unqframes(k) - (memory-1)));
    % so if curframes is 15 15 13 14 15 15 15, the present frame is 15, and
    % memory==1, cullind will return [3 4].  If memory==2, returns [3].
    
    % DO NOT remove lost objects from the objs matrix -- keep them, with a
    % track ID that simply ends at this frame.
    % objs(6, cur(cullind)) = 0; % particles with trackid 0 are culled at end
    cur(cullind) = [];
    
    % add the rest to leftovers -- the set of indices of objects not linked
    % in recent earlier frames
    leftovers = cur;
    
    if diagnostic
        if ~isempty(leftovers)
            disp('Diagnostic: IDs of Leftovers')
            objs(6,leftovers)
            curframes
            (unqframes(k) - (memory-1))
            cullind
        end
    end
    
    % assign new trackids to members of nex for which bmap ~= inv(fmap)
    % complicated -- from Andy
    bmapind = (1:length(bmap))';
    onesnew = ones(size(bmapind));
    onesinv = ones(size(notlost));
    resinv = kron(fmap(notlost), onesnew);
    resnew = kron(onesinv, bmapind);
    mask = all(resinv ~= resnew, 2);
    newguys = bmapind(mask);

    if (~isempty(newguys))
        objs(6, nex((newguys))) = nextid:(nextid+length(newguys)-1);
        nextid = nextid + length(newguys);
    end
    
    % show progress
    if dispopt
        waitbar(k/(length(unqframes)-1), progbar, progtitle);
    end
end
if dispopt
    close(progbar)
end

% now get rid of unmatched particles
singletons = find(objs(6,:) == 0);
objs(:, singletons) = [];

if dispopt
    % Display information about tracks found
    ntracks = length(unique(objs(6,:)));
    fs = sprintf('Found %d unique tracks.', ntracks);
    disp(fs)
    if ntracks <= 10
        for j=unique(objs(6,:))
            objs_j = objs(:,objs(6,:)==j);  % track j
            frj = objs_j(5,:);
            x = objs_j(1,:);  % x positions of this track
            y = objs_j(2,:);  % y positions of this track
            stdj = sqrt(var(x) + var(y));  % standard deviation of this track's position
            %frj = objs(5,find(objs(6,:)==j));
            fs = sprintf('   Track %d: frames %d:%d (%d frames), std. %.1f pixels',...
                j, min(frj), max(frj), length(frj), stdj);
            % Note that some frames between min and max may not be good
            disp(fs)
        end
    else
        disp('More than 10 tracks found; not displaying details.');
    end
end

if memory>1
    disp('nnlink_rp:')
    disp('Warning: "memory>1" works properly only for continuous frame numbers.')
    disp('   Else, adjacent frames are considered even if frame numbers are')
    disp('   not contiguous.  Try, e.g, with frames [3 4 5 11 12 13] -- ')
    disp('   5 and 11 will be cur and nex frames, though they''re separated');
end


