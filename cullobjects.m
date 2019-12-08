% cullobjects.m
%
% function to cull objects from an object array based on their width and rms
% distance of gradient lines to center (each calculated in fo5_rp.m)
%
% Note that the culling is applied to individual objects in
% individual frames, and so does not consider "track" information.
% It may, therefore, generate tracks with objects missing in frames.
%
% Inputs
% objs : object array, linked or not, e.g. output by fo5_rp.m
%        Note that row 7 is the object width (sigma), and row 8 is the mean
%            distance-*squared* between gradient lines and center for
%            radial-symmetry localization.  (If another localization method
%            was used, the row 8 values will all be zero and no culling of
%            these will be done.)
% cutoffs = [sigmamaxstd dmaxstd] cutoff values, # std. above median
%         = [] use automatic cutoffs for sigma and d (each median + 0.3 std. dev.)
%         = -1 : ask for user input (if cutoffs is any negative number)
% culloption : culling method
%            'rect' : [default] based on sigma > threshold and d >
%               threshold, separately, like rectangles in the sigma, d space
%              'diag' : cull based on 
%               d > (0.5*(sigmamaxcutoff - sigma) + dmaxcutoff),
%               like a diagonal cut in the sigma, d space.  See notes Feb
%               2013.  Only allowed if row 8 has nonzero values.
%
% plotopt : if true, plot the histograms of sigma and sqrt(d2min)
% doasses : call assess_localization.m to assess localization accuracy, if
%           ground-truth values are known.
% fluos : ground truth array (leave empty to ignore)
%
% Output
% objs_cut : object array in which the culled columns have been deleted.
% sigmamaxstd : sigma cutoff, # std. above median
% dmaxstd : sigma cutoff, # std. above median
%         
% Written originally for the 2013 "Localization Challenge"
%
% Raghuveer Parthasarathy
% Feb. 27, 2013
% last modified April 18, 2013


function [objs_cut sigmamaxstd dmaxstd] = cullobjects(objs, cutoffs, culloption, plotopt, doassess, fluos)

if ~exist('plotopt', 'var') || isempty(plotopt)
    plotopt = false;
end
% But: force plotting if user-input is requested
if min(cutoffs) < 0
    plotopt = true;
end

if ~exist('doassess', 'var') || isempty(doassess)
    doassess = false;
end


% Culling option; default is radial symmetry (RP 2012)
if ~exist('culloption', 'var') || isempty(culloption)
    culloption = 'rect';
end


objs_cut = objs;

if doassess
   assess_localization(objs_cut, fluos);
end


%% width, sigma
sigma = objs_cut(7,:);
[ns, binns] = hist(sigma,200);
means = mean(sigma);
meds = median(sigma);
stds = std(sigma);
peaks = mean(binns(ns==max(ns(:))));  % mean in case there's more than one equal max.

if plotopt
    hsigma = figure; plot(binns, ns, 'x-', 'color', [0.9 0.6 0.2])
    xlabel('Width (Sigma), px');
    hold on
end
fs = sprintf('   Width:  peak %.3f, mean %.3f, median %.3f, std %.3f', ...
    peaks, means, meds, stds); disp(fs)
if isempty(cutoffs)
    sigmamaxcutoff = meds + 0.3*stds;
    fs = sprintf('   Using max. width cutoff %.3f px', sigmamaxcutoff); disp(fs);
elseif min(cutoffs) < 0
    fs = sprintf('   Median + 0.3 std. = %.3f px', meds + 0.3*stds); disp(fs)
    sigmamaxcutoff = input('Max cutoff (absolute, not #std!): ');  % (peakbr+2*stdbr)
else
    sigmamaxcutoff = meds + cutoffs(1)*stds;
end
sigmamaxstd = (sigmamaxcutoff-meds)/stds;  % for output

%% Goodness of radial-symmetry-fit (mean of sqrt(distance^2))
d = sqrt(objs_cut(8,:));
if max(d) > 0
    % there is d information (see comments header)
    cull_d = true;
    [nd, binnd] = hist(d,200);
    meand = mean(d);
    medd = median(d);
    stdd = std(d);
    peakd = mean(binnd(nd==max(nd(:))));  % mean in case there's more than one equal max.
    if plotopt
        hd = figure; plot(binnd, nd, 'x-', 'color', [0.2 0.6 0.8])
        hold on
        xlabel('Gradient line min distance, px');
    end
    
    fs = sprintf('   d to center:  peak %.3f, mean %.3f, median %.3f, std %.3f', ...
        peakd, meand, medd, stdd); disp(fs)
    
    if isempty(cutoffs)
        dmaxcutoff = medd + 0.3*stdd;
        fs = sprintf('    Using max. d-to-ctr cutoff %.3f px', dmaxcutoff); disp(fs);
    elseif min(cutoffs) < 0
        fs = sprintf('   Median + 0.3 std. = %.3f px', medd + 0.3*stdd); disp(fs)
        dmaxcutoff = input('Max cutoff (absolute, not #std!): ');  % (peakbr+2*stdbr)
    else
        dmaxcutoff = medd + cutoffs(2)*stdd;
    end
    dmaxstd = (dmaxcutoff-medd)/stdd;  % for output
    switch lower(culloption)
        case {'rect'}
            % Cull (rectangular cuts)
            badsd = sigma>sigmamaxcutoff | d>dmaxcutoff;
        case {'diag'}
            % Cull (diagonal cut)
            badsd = d > (0.5*(sigmamaxcutoff - sigma) + dmaxcutoff);
    end
else
    % there is d information (see comments header)
    cull_d = false;
    dmaxcutoff = 9e99;
    dmaxstd = 9e99;
    % Cull (sigma only)
    badsd = sigma>sigmamaxcutoff;
end    



% Cull!
objs_cut(:,badsd)=[];

if plotopt
    figure(hsigma)
    plot([sigmamaxcutoff sigmamaxcutoff], [0 max(ns)], 'k:')
    if cull_d
        figure(hd)
        plot([dmaxcutoff dmaxcutoff], [0 max(nd)], 'k:')
    end
end

if doassess
   assess_localization(objs_cut, fluos);
end

