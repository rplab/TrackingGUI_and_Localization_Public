% TrackingGUI_rp.m
% 
% Copyright 2012, Raghuveer Parthasarathy, The University of Oregon
%
%%
% Disclaimer / License  
%   This program is free software: you can redistribute it and/or 
%     modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.
%   This set of programs is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
%   General Public License for more details.
%   You should have received a copy of the GNU General Public License 
%   (gpl.txt) along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% TrackingGUI_rp.m: 
% A graphical interface for particle tracking, allowing 
%   -- choice of tracking algorithm
%   -- linkage of objects -> tracks
%   -- display of tracking output overlayed on image frames
%
% Inputs: 
%    im : (optional) 3D array of 2D images, e.g. from TIFFseries.m
%
% Outputs:
%    *None* (!).  It's very difficult to get output from a GUI.  (See various web
%    pages.  Will just save arrays and parameters to a MAT file:
%       objsize : object size parameter used for particle localization
%       thresh : threshold parameter used for particle localization
%       objs : object matrix, output by im2obj_rp.m (see file)
%       objs_link : linked object matrix, output by nnlink_rp.m (see file)
%    Or: save 'simple' output as a text file -- just x and y positions of
%       each linked object, without additional information.
%       Row 1 = x positions (px) of object #1 in each frame in which it
%       exists; Row 2 = y positions.
%       Row 3 = x positions of object #2; Row 4 = y positions. 
%       etc.
%       Note that the columns correspond to frames, but the set of
%       frames in which each object exists need not be the same. 
%       The tab-delimited text file can be opened e.g. by WordPad or Excel
%
% TrackingGUI_rp.m  begun April 3, 2012.  
%
% The display and transparency functions as well as the object information
% selection / data table are taken from WaterGUI[1-6].m by Matt Jemielita, 2010
%
% Raghuveer Parthasarathy
% Modifications since June 2012:
% June 10, 2012 (added text  output option)
% June 27, 2012 (fixed bug for frame slider range if single image is input;
%                consistent definition of threshold default values)
% August 25, 2012: uses fo5_rp.m and allows input of prior neighborhood
%    positions, gradient voting for region finding, etc.
% September 9, 2016: reconcile with Tristan Hormel's version, including
%    Tristan's April 2014 routines for visualizing results of phase 
%    separated domain tracking (bilateral filter, etc.) Note that watershed
%    segmentation code is commented out (hdispwater, etc.)
% Dec. 17, 2018: adjust GUI size
% July 17, 2019: allow adjustable display intensity (sliders)
% Last modified Sept. 27, 2019: minor edit in line 1205, to allow empty
%    tmpobj (no objects)

function fGUI = TrackingGUI_rp(im)

% A nested function

%% Get Filenames and info for loading images
% -----------------------------------------------------------------------
programdir = cd; % present directory.

if ~exist('im', 'var') || isempty(im)
    % images were not input; get file name info
    imagesloaded = false;
    %    Can be multipage TIFF
    [fbase, frmin, frmax, formatstr, FileName1, ~, PathName1, ext, ismultipage] = ...
        getnumfilelist;
    % If there's only one image, getnumfilelist returns empty arrays for
    % frmin, frmax, formatstr.  Replace frmin and frmax with 1s.
    if isempty(frmin)
        frmin = 1;
        frmax = 1;
    end
else
    % images were input
    imagesloaded = true;
    frmin = 1;
    frmax = size(im,3);
    FileName1 = '[input 3D array]';
end
Nframes = frmax - frmin + 1;

% -----------------------------------------------------------------------
%% GUI components
% -----------------------------------------------------------------------

%  Initialize and hide the GUI as it is being constructed.
scrsz = get(0,'ScreenSize');  % screen size

fGUI = figure('Name','TrackingGUI_rp', 'Menubar','none', ...
    'Visible','off','Position',[1 1 scrsz(3) scrsz(4)], ...
    'Color', [0 40 90]/255);  
            
% Construct the components.  Create axes, for the images
% Initialization
% Initialize file uploading controls. All these controls will be contained
% in the hupdate panel.
% hfileupdate controls what files will be analyzed by the program

% Signature
uicontrol('Style','text','BackgroundColor', [0.5 0.7 0.7], ...
    'String','TrackingGUI_rp.m.  Raghuveer Parthasarathy, 2012-', ...
    'FontWeight', 'bold', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',...
    [0.67 0.96 0.25 0.03]);
%Create an exit button
hexit = uicontrol('Style','pushbutton',...
    'String','Exit','Units', 'normalized','FontWeight', 'bold',...
    'Position',[0.94 0.86 0.05 0.1], 'Callback',{@exit_Callback});

%% Panel for loading the files that will be analyzed by the GUI 
himageselect = uipanel('Title', 'Display Frame', 'FontSize', 11, ...
    'Units', 'normalized', 'Position', [0.67 0.86 0.25 0.1]);
   % Primary frame no. to load, analyze
% hFileNameText = 
uicontrol('Parent', himageselect, ...
    'Style', 'text', 'String', FileName1,...
    'FontAngle', 'italic', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.55 0.9 0.4]);
% hframenotextinit = 
uicontrol('Parent', himageselect,'Style','text','Units',...
    'normalized','String','Frame no.', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.1 0.2 0.4]);
hframenotext  = uicontrol('Parent', himageselect,'Style','edit','Units', ...
    'normalized', 'Position',[0.35, 0.1, 0.15, 0.4], ...
    'Callback',{@framenotext_Callback});
hframeno = uicontrol('Parent', himageselect,'Style','slider', ...
    'Max', frmax + 0.1*(Nframes==1), 'Min', frmin, 'Value', frmin, 'Units', 'normalized',...
    'SliderStep', [1/Nframes min([0.1, 10/Nframes])], 'Position', ...
    [0.55 0.1 0.35 0.4], 'Callback',{@frameno_Callback});

%% Panel for adjusting display range
hAdjustDisplay = uipanel('Title', 'Level', 'FontSize', 11, ...
    'Units', 'normalized', 'Position', [0.55 0.70 0.05 0.29]);
% Automatic display intensity range
hAutoDisplay = uicontrol('Parent', hAdjustDisplay,'Style','checkbox',...
    'String','Auto','Fontweight', 'bold', 'Units', 'normalized', ...
    'Position',[0.05, 0.8,0.9,0.15],'Value', 1, 'Callback',{@displayRange_Callback});
% Display intensity range (0-1)
hDisplayRangeMin = uicontrol('Parent', hAdjustDisplay,'Style','slider', ...
    'Max', 1, 'Min', 0, 'Value', 0, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', ...
    [0.05 0.05 0.4 0.7], 'Callback',{@displayRange_Callback});
hDisplayRangeMax = uicontrol('Parent', hAdjustDisplay,'Style','slider', ...
    'Max', 1, 'Min', 0, 'Value', 1, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', ...
    [0.55 0.05 0.4 0.7], 'Callback',{@displayRange_Callback});


%% Neighborhood finding parameters
hnhoodparampanel = uipanel('Title','Neighborhood Parameters','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .57 .32 .25]);
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','Process Option', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.02 0.87 0.27 0.12]);
% processing options for fo5_rp.m -- note that these strings need to be
% identical to those expected by fo5_rp.m
processoptarray = {'spatialfilter', 'gradientvote', 'none'};
hprocessopt = uicontrol('Parent', hnhoodparampanel, 'Style','popupmenu', ...
    'String', strcat(char(processoptarray(1)), ' | ', char(processoptarray(2)), ' | ', ...
    char(processoptarray(3))), 'Units', 'normalized','Position', [0.3,0.87,0.3,0.12], ...
    'Callback',{@processopt_Callback});

% Spatial filtering options
defaultobjsize = 7;
spatialfilt_ypos = 0.71;  % (relative) y-position of spatial filtering option parameters
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','bpfiltsize', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.02 spatialfilt_ypos 0.15 0.1]);
hbpfiltsizetext  = uicontrol('Parent', hnhoodparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%d', defaultobjsize), ...
    'Position',[0.17, spatialfilt_ypos, 0.08, 0.1], ...
    'Callback',{@bpfiltsize_Callback}); 
hbpfiltsize = uicontrol('Parent', hnhoodparampanel,'Style','slider', ...
    'Max', 100, 'Min', 0, 'Value', defaultobjsize, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', [0.26 spatialfilt_ypos 0.19 0.1], ...
    'Callback',{@bpfiltsize_Callback}); 
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','nsize', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.47 spatialfilt_ypos 0.1 0.1]);
hnsizetext  = uicontrol('Parent', hnhoodparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%d', defaultobjsize), ...
    'Position',[0.58, spatialfilt_ypos, 0.08, 0.1], ...
    'Callback',{@nsize_Callback});
hnsize = uicontrol('Parent', hnhoodparampanel,'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultobjsize, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', [0.67 spatialfilt_ypos 0.19 0.1], ...
    'Callback',{@nsize_Callback});
% Button to lock filtering object size and neighborood object size together
hlockobjsize  = uicontrol('Parent', hnhoodparampanel,'Style','toggle',...
    'String','Lock','Units', 'normalized', 'BackgroundColor', [0.8 1.0 0.6], ...
    'Position',[0.89, spatialfilt_ypos,0.1,0.1],'Value', 1, ...
    'Callback',{@lockobjsize_Callback});

% Gradient voting options
defaultobjsize = 7;
gradvote_ypos = 0.57;  % (relative) y-position of gradient voting option parameters
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','gradobjsize', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.02 gradvote_ypos 0.15 0.1]);
hgrobjsizetext  = uicontrol('Parent', hnhoodparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%d', defaultobjsize), 'Position',[0.17, gradvote_ypos, 0.08, 0.1], ...
    'Callback',{@grobjsize_Callback});
hgrobjsize = uicontrol('Parent', hnhoodparampanel,'Style','slider', ...
    'Max', 100, 'Min', 0, 'Value', defaultobjsize, 'Units', 'normalized',...
    'SliderStep', [0.01 0.1], 'Position', [0.26 gradvote_ypos 0.19 0.1], ...
    'Callback',{@grobjsize_Callback});
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','dir.', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.47 gradvote_ypos 0.05 0.1]);
hgraddiropt = uicontrol('Parent', hnhoodparampanel, 'Style','popupmenu', ...
    'String', 'Both dir.|Positive|Negative', ...
    'Units', 'normalized','Position', [0.53,gradvote_ypos,0.21,0.12], ...
    'Callback',{@setprocessparam});
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','grthresh', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.76 gradvote_ypos 0.12 0.1]);
defaultgrthresh = 0.99;
hgrthreshtext  = uicontrol('Parent', hnhoodparampanel,'Style','edit','Units', ...
    'normalized', 'String', sprintf('%.3f', defaultgrthresh), ...
    'Position',[0.88, gradvote_ypos, 0.11, 0.08], ...
    'Callback',{@setprocessparam});

% hthreshoptinit = 
uicontrol('Parent', hnhoodparampanel,'Style','text','Units',...
    'normalized','String','Threshold Option', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.29 0.3 0.12]);
hthreshopt = uicontrol('Parent', hnhoodparampanel, 'Style','popupmenu', ...
    'String', '1. Intens. Thresh.|2. Std. Thresh.|3. Bright N', ...
    'Units', 'normalized','Position', [0.05,0.16,0.3,0.12], ...
    'Callback',{@threshopt_Callback});

%hthreshtextinit1 = 
uicontrol('Parent', hnhoodparampanel, 'Style','text',...
    'String','thr (0-1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.32,0.13,0.1]);
defaultthresh1 = 0.99;
hthresh1 = uicontrol('Parent', hnhoodparampanel, 'Style','slider', ...
    'Max', 0.99999, 'Min', 0.00, ...
    'SliderStep', [0.02 0.1], 'Units', 'normalized',...
    'Value', defaultthresh1, 'Position', [0.7,0.32,0.27,0.1], ...
    'Callback',{@thresh1_Callback});
hthreshtext1  = uicontrol('Parent', hnhoodparampanel, 'Style','edit',...
    'Units', 'normalized', 'String', sprintf('%.4f', get(hthresh1, 'Value')), ...
    'Position',[0.53,0.32,0.15,0.1], ...
    'Callback',{@threshtext1_Callback});
%hthreshtextinit2 = 
uicontrol('Parent', hnhoodparampanel, 'Style','text',...
    'String','thr (>1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.17,0.13,0.1]);
defaultthresh2 = 3.0;
hthreshtext2  = uicontrol('Parent', hnhoodparampanel, 'Style','edit',...
    'String', sprintf('%.1f', defaultthresh2), 'Units', 'normalized', ...
    'Position',[0.53,0.17,0.15,0.1], ...
    'Callback',{@threshtext2_Callback});
hthresh2 = uicontrol('Parent', hnhoodparampanel, 'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultthresh2, ...
    'SliderStep', [0.01 0.10], 'Units', 'normalized',...
    'Position', [0.7,0.17,0.27,0.1], ...
    'Callback',{@thresh2_Callback});
defaultthresh3 = 3;
%hthreshtextinit3 = 
uicontrol('Parent', hnhoodparampanel, 'Style','text',...
    'String','N (>=1)', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.38,0.02,0.13,0.1]);
hthreshtext3  = uicontrol('Parent', hnhoodparampanel, 'Style','edit',...
    'String', sprintf('%d', defaultthresh3), 'Units', 'normalized', ...
    'Position',[0.53,0.02,0.15,0.1], ...
    'Callback',{@threshtext3_Callback});
hthresh3 = uicontrol('Parent', hnhoodparampanel, 'Style','slider', ...
    'Max', 100, 'Min', 1, 'Value', defaultthresh3, ...
    'SliderStep', [0.01 0.10], 'Units', 'normalized',...
    'Position', [0.7,0.02,0.27,0.1], ...
    'Callback',{@thresh3_Callback});

%% Displaying filtered and thresholded images
hdispprocesspanel = uipanel('Title','Display Processed Imgs','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .46 .32 .1]);
% hdispimgText = 
uicontrol('Parent', hdispprocesspanel, ...
    'Style', 'text', 'String', 'Processed for finding maxima, not for fitting',...
    'FontAngle', 'italic', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position',[0.05 0.7 0.9 0.3]);
hdispprocess  = uicontrol('Parent', hdispprocesspanel,'Style','checkbox',...
    'String','DispProcess','Units', 'normalized', ...
    'Position',[0.05, 0.2,0.3,0.4],'Value', 0,...
    'Callback',{@dispprocess_Callback});
hdispthreshprocess    = uicontrol('Parent', hdispprocesspanel,'Style','checkbox',...
    'String','DispThreshProcess','Units', 'normalized', ...
    'Position',[0.4, 0.2,0.35,0.4],'Value', 0,...
    'Callback',{@dispthreshprocess_Callback});
% hdispwater = uicontrol('Parent',hdispprocesspanel,'Style','checkbox',...
%     'String','DispWatershed','Units', 'normalized', ...
%     'Position',[0.8, 0.2, 0.2, 0.4],'Value',0,...
%     'Callback',{@dispwater_Callback});

%% Tracking
htrackpanel = uipanel('Title','Tracking','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .18 .32 .27]);

% Center-finding algorithm option
ctrfinding_ypos = 0.8;
uicontrol('Parent', htrackpanel, 'Style','text',...
    'String','Localiz. method', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.02,ctrfinding_ypos,0.23,0.12]);
ctrfitstrarray = {'radial', 'gaussmle', 'nonlineargauss', 'lineargauss','centroid'};
hctrfitstr = uicontrol('Parent', htrackpanel, 'Style','popupmenu', ...
    'String', strcat(char(ctrfitstrarray(1)), ' | ', ...
                     char(ctrfitstrarray(2)), ' | ', ...
                     char(ctrfitstrarray(3)), ' | ', ...
                     char(ctrfitstrarray(4)), ' | ', ...
                     char(ctrfitstrarray(5))), ...
    'Units', 'normalized','Position', [0.27,ctrfinding_ypos+0.01,0.2,0.12],...
    'Value', 1, 'Callback',{@ctrfitstr_Callback});
husenhoodctrs = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Use prev. ctrs.','Fontweight', 'bold', 'Units', 'normalized', ...
    'Position',[0.49, ctrfinding_ypos+0.01,0.28,0.12],'Value', 0);
h1pernhood = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','1/nhood','Fontweight', 'bold', 'Units', 'normalized', ...
    'Position',[0.79, ctrfinding_ypos,0.19,0.12],'Value', 0);

% Orientation-finding algorithm option
uicontrol('Parent', htrackpanel, 'Style','text',...
    'String','Orient. method', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.02,0.65,0.23,0.12]);
orientstrarray = {'none', 'momentcalc'};
horientstr = uicontrol('Parent', htrackpanel, 'Style','popupmenu', ...
    'String', strcat(char(orientstrarray(1)), ' | ',char(orientstrarray(2))), ...
    'Units', 'normalized','Position', [0.27,0.65,0.2,0.12],...
    'Value', 1, 'Callback',{@orientstr_Callback});

% Size-finding algorithm option (Tristan Hormel; for vesicle domains)
uicontrol('Parent', htrackpanel, 'Style','text',...
    'String','Size method', 'FontWeight', 'bold', 'Units', 'normalized',...
    'Position',[0.5,0.65,0.2,0.12]);
sizestrarray = {'none','bilateral_filter'};
hsizestr = uicontrol('Parent',htrackpanel, 'Style','popupmenu', ...
    'String', strcat(char(sizestrarray(1)), ' | ', ...
                     char(sizestrarray(2))), ...
    'Units', 'normalized','Position', [0.7 0.65 0.2 0.12],...
    'Value', 1, 'Callback',{@sizestr_Callback});

% Tracking
htrackthisfr    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Track this frame','Units', 'normalized',...
    'Position',[0.05,0.47,0.3,0.12],'value', 0, ...
    'Callback',{@trackthisfr_Callback});
htrackthisdone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.47,0.3,0.12],'value', 0, ...
    'Callback',{@trackthisdone_Callback});
htrackallfr    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Track all frames','Units', 'normalized',...
    'Position',[0.05,0.32,0.3,0.12],'value', 0, ...
    'Callback',{@trackallfr_Callback});
htrackalldone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.32,0.3,0.12],'value', 0, ...
    'Callback',{@trackalldone_Callback});
hcleartracking    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Clear','Units', 'normalized',...
    'Position',[0.55,0.32,0.15,0.12],'value', 0, ...
    'Callback',{@cleartracking_Callback});

linkopt_ypos = 0.15;
uicontrol('Parent', htrackpanel,'Style','text','Units',...
    'normalized','String','Link MaxStep^2', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.05 linkopt_ypos 0.25 0.12]);
hlinksteptext  = uicontrol('Parent', htrackpanel,'Style','edit','Units', ...
    'normalized', 'Position',[0.32, linkopt_ypos, 0.1, 0.12], 'String', '1000', ...
    'Callback',{@linksteptext_Callback});
uicontrol('Parent', htrackpanel,'Style','text','Units',...
    'normalized','String','Link Memory', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[0.45 linkopt_ypos 0.2 0.12]);
hlinkmemorytext  = uicontrol('Parent', htrackpanel,'Style','edit','Units', ...
    'normalized', 'Position',[0.67, linkopt_ypos, 0.1, 0.12], 'String', '0', ...
    'Callback',{@linkmemorytext_Callback});
hlink    = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Link Tracks','Units', 'normalized',...
    'Position',[0.05,0.02,0.3,0.1],'value', 0, ...
    'Callback',{@link_Callback});
hlinkdone    = uicontrol('Parent', htrackpanel,'Style','checkbox',...
    'String','Done','Units', 'normalized',...
    'Position',[0.4,0.02,0.3,0.1],'value', 0, ...
    'Callback',{@linkdone_Callback});
hclearlink = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','Clear','Units', 'normalized',...
    'Position',[0.55,0.02,0.15,0.1],'value', 0, ...
    'Callback',{@clearlink_Callback});

hculllink = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','CullObjs','Units', 'normalized',...
    'Position',[0.8,0.15,0.15,0.15],'value', 0, ...
    'Callback',{@culllink_Callback});
hunculllink = uicontrol('Parent', htrackpanel,'Style','pushbutton',...
    'String','un-Cull','Units', 'normalized',...
    'Position',[0.8,0.02,0.15,0.1],'value', 0, ...
    'Callback',{@unculllink_Callback});

%% Save output, in a MAT file, or load previously calculated values
hsavepanel = uipanel('Title','Save / Load','FontSize',11,...
'Units',  'normalized', 'Position',[.67 .05 .32 .10]);
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Save output','Units', 'normalized',...
    'Position',[0.05,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@saveoutput_Callback});
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Save TXT','Units', 'normalized',...
    'Position',[0.35,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@savetxt_Callback});
uicontrol('Parent', hsavepanel,'Style','pushbutton',...
    'String','Load obj data','Units', 'normalized',...
    'Position',[0.65,0.1,0.25,0.8],'value', 0, ...
    'Callback',{@loadoutput_Callback});

%% Messages
hmsg    = uicontrol('Style','text','BackgroundColor', [1 1 0.7],...
    'String','Message: ', 'FontWeight', 'bold', 'Units', 'normalized',...
    'HorizontalAlignment', 'left', 'Position', [0.5 0.05 0.15 0.12]);

%% Display tracks, etc.

hdisptrackspanel = uipanel('Title','Display Tracks','FontSize',11,...
'Units',  'normalized', 'Position',[.35 .05 .11 .15]);
hdispcircles    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show objs','Units', 'normalized', ...
    'Position',[0.05, 0.6,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});
hdispIDs    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show IDs','Units', 'normalized', ...
    'Position',[0.05, 0.3,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});
hdisptracks    = uicontrol('Parent', hdisptrackspanel,'Style','checkbox',...
    'String','Show Tracks','Units', 'normalized', ...
    'Position',[0.05, 0.01,0.8,0.25],'Value', 0,...
    'Callback',{@displayimage});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the button group which allows the user to either collect
% information about particular cells in the image or to have a zoom feature
% for the image shown.
hbgroup = uibuttongroup('Title', 'Control Display ','visible','off',...
    'Position',[0.23 0.05 .1 .15]);
% Create three radio buttons in the button group.  (Could name the handles
% to the controls, but not necessary)
uicontrol('Style','Radio','String','Zoom','Units', 'normalized',...
    'pos',[.10 .05 .90 .25],'Parent',hbgroup, ...
    'HandleVisibility','off');
uicontrol('Style','Radio','String','Pan','Units',...
    'normalized', 'pos',[.10 .35 .900 .25],'Parent',hbgroup,...
    'HandleVisibility','off');
uicontrol('Style','Radio','String','Object Information','Units',...
    'normalized', 'pos',[.10 .65 .900 .25],'Parent',hbgroup,...
    'HandleVisibility','off');
% Initialize some button group properties. 
set(hbgroup,'SelectionChangeFcn',@selcbk);
set(hbgroup,'SelectedObject',[]);  % No selection
set(hbgroup,'Visible','on');
%Initializing the datacursormode
dcm_obj = datacursormode(fGUI);
%Initially the datacursormode() will be set to off. By clicking on the
%object information button it will be turned on.
set(dcm_obj,'Enable','off');
set(dcm_obj, 'Updatefcn', @ObjectInfo_Callback);

% ------------------------------------------------------------------
%% Initialize the GUI.

% ------------------------------------------------------------------
% Create all variables here

% Outputs
% NOT output in the function (very hard to do this in a GUI).  After 
% calculation, click 'save' to write to MAT file
processopt = char(processoptarray(get(hprocessopt,'Value')));
processparam = [];
lockobjsize = get(hlockobjsize, 'Value');
thresh = get(hthresh1, 'Value');
threshopt = [];
% Call setprocessparam, so that present (default) parameter values are
% in-place
setprocessparam;

objs = [];
objs_link = [];
savedobjs = [];   % to revert if un-culling objects
savedobjs_link = []; % to revert if un-culling objects

% Other
ctrfitstr = char(ctrfitstrarray(get(hctrfitstr,'Value')));
orientstr = char(orientstrarray(get(horientstr,'Value')));
sizestr = char(sizestrarray(get(hsizestr,'Value')));
fitstr = {ctrfitstr; orientstr; sizestr};
lsqoptions = optimset('lsqnonlin');
linkstep = round(str2double(get(hlinksteptext,'string')));
linkmem =  round(str2double(get(hlinkmemorytext,'string')));

% Load one image, or use the first frame of the 'im' array
% In general, will use 'A' for the displayed image
if imagesloaded
    A = im(:,:,1);
    isColor = false;
else
    cd(PathName1) %Go to the directory of the images
    if ~ismultipage
        A=imread(FileName1);    % load the first file image
    else
        % multipage TIFF
        A = imread(FileName1, 1);
    end
    isColor = (ndims(A)==3);  % TRUE if the image is color (3 layers)
    if isColor
        prompt = {'If color: which channel to use?  (1, 2, 3): '};
        dlg_title = 'Color option'; num_lines= 1;
        % guess at which channel to use as default, based on total brightness
        b = zeros(1,3);
        for ch = 1:3
            b(ch) = sum(sum(A(:,:,ch)));
        end
        [~, ic] = max(b);
        def     = {num2str(ic)};  % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        colchannel = round(str2double(answer(1)));
        A = A(:,:,colchannel);  % select the color channel to examine
        pause(0.5);  % seems to be necessary to keep MATLAB from crashing!
    end
end
% Determine the bit depth of the image, for adjustable display range
if isa(A, 'uint8') % is A a 'uint8' variable?
    bitDepth = 8;
elseif isa(A, 'uint16') % is z a 'uint16' variable?
    bitDepth = 16;
else
    bitDepth = [];
end
displayRange = []; % 0-1, Define here to use later. Max caxis will be displayRange*2^bitDepth-1)

currframe = frmin;  % # of the current (primary) frame

% Show filtered, thresholded image
showprocess = false;
showthreshprocess = false;
% showwatershed = false;
showbilateralfilter = false;
processedA = [];  % filtered image
threshprocessA = [];  % thresholded, filtered image
markedA = [];
bilateralA = [];

cd(programdir); % Go back to the directory from which the GUI was called

% Handles to display images
% Am I going to use these?
im1 = [];  % handle to primary image display

% Status variables
trackthisdone = false(Nframes,1);  % true for frames that have been segmented
islinkdone = false;  % is linkage of objects into tracks done?

% necessary?
colormap('gray');
        
outfile = [];  % file name for saving results (MAT file)

% Creates axes on the figure. At a later time in the program these axes will
% be filled with the images we are analyzing.
imageRelHeight = 0.78;
imageRelWidth = size(A, 2)/size(A,1)*imageRelHeight;
axeshandle = axes('Position', [0.01 0.22 imageRelWidth imageRelHeight]);

%Creating a table where the information about the cell we click on is
%displayed.
selectedobjectdata = cell(6,1);
t = uitable('Parent', fGUI, 'Units', 'normalized','Position',...
   [0.01 0.05 0.2 0.15]);
set(t, 'ColumnName', []);
set(t, 'RowName',{ 'Centroid x',  'Centroid y', 'Brightness',...
    'Particle ID', 'Frame no.', 'Track ID', 'Sigma', 'mean d^2'});
set(t, 'Data', selectedobjectdata);
align([t axeshandle], 'Fixed',100,  'Top');

% Move the GUI to the center of the screen.
movegui(fGUI,'center')
% Make the GUI visible.
set(fGUI,'Visible','on');
set([htrackthisfr, htrackallfr, hlink, hculllink, hunculllink, ...
    hexit, hmsg], 'BackgroundColor',[0.85 1.0 0.6]);
set([hcleartracking, hclearlink], 'BackgroundColor',[1.0 0.7 0.4]);
set([hexit, hmsg], 'BackgroundColor',[1.0 0.5 0.2]);

figure(fGUI);

% Set Frame number slider and text.
set(hframenotext, 'String', sprintf('%d', currframe));
set(hframeno, 'Value', currframe);
% Call the processing and threshold option functions, just to highlight the
% default rows
processopt_Callback;
threshopt_Callback;

% Display the first image
displayimage

%% ----------------------------------------------------------------------
% **********************************************************************

% Callbacks and other functions

    function exit_Callback(source,eventdata)
        % Exit
        close all
    end

% Callback functions for loading and displaying frames

    function framenotext_Callback(source,eventdata)
        % set frame to view (text entry)
        currframe = round(str2double(get(source,'string')));
        % Also update slider:
        set(hframeno, 'Value', currframe);
        % Load the image
        loadimages;
        % Display
        displayimage;
        % Update the checkbox indicating whether segmentation has been done
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function frameno_Callback(hObject, source,eventdata)
        % set primary frame to view (slider)
        currframe = round(get(hObject,'Value'));
        % Also update text entry box:
        set(hframenotext, 'String', sprintf('%d', currframe));
        % Load the image
        loadimages;
        % Display
        displayimage;
        % Update the checkbox indicating whether segmentation has been done
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function loadimages(source,eventdata)
        set(hmsg, 'String', 'Loading images...');
        if imagesloaded
            % all images are in the 'im' array
            A = im(:,:,currframe);
        else
            % load from file
            cd(PathName1)     %Go to the directory of the images
            % load the image into "A"
            % Don't do any calculations or updates
            % Note that all file names were determined previously
            if ~ismultipage
                framestr = sprintf(formatstr, currframe);
                A  = imread(strcat(fbase, framestr, ext));
            else
                % multipage TIFF
                A  = imread(strcat(fbase, ext), currframe);
            end
            cd(programdir) %go back to the original directory
            if isColor
                % select the color channel to examine
                A = A(:,:,colchannel);
            end
        end
        set(hmsg, 'String', 'Image loaded.');
        % if the boxes are checked, calculate the filtered and thresholded
        % images
        if showprocess || showthreshprocess
            processedA = calcprocessimg;
        end
        if showthreshprocess
            threshprocessA = calcthreshprocess;
        end
%         if showwatershed
%             markedA = calcmarkers;
%         end
        if showbilateralfilter
            bilateralA = calcbilateralfilter;
        end
    end

% Set Display range
    function displayRange_Callback(hObject, source,eventdata)
        if ~get(hAutoDisplay, 'Value')
            % Unchecked, so use slider values
            % Check that min < max; adjust
            if get(hDisplayRangeMin, 'Value') >= get(hDisplayRangeMax, 'Value')
                set(hDisplayRangeMin, 'Value', 0.99*get(hDisplayRangeMax, 'Value'))
            end
            if get(hDisplayRangeMax, 'Value') <= get(hDisplayRangeMin, 'Value')
                set(hDisplayRangeMax, 'Value', 1.01*get(hDisplayRangeMin, 'Value'))
            end
            displayRange = [get(hDisplayRangeMin, 'Value') get(hDisplayRangeMax, 'Value')];
        else
            % Checked, so make empty, to AutoRange
            displayRange = [];
        end
        displayimage
    end

% Processing for neighborhood determination parameters (Object size, etc.)
% and threshold parameters

    function processopt_Callback(source,eventdata)
        % set processing option, and highlight accordingly
        % also call setprocessparam to set processing parameter values
        ltgray = 0.95*[1 1 1];
        highlightcolor = [1 1 0.7];
        colorarray = [ltgray; highlightcolor];
        processopt = char(processoptarray(get(hprocessopt,'Value')));
        switch processopt
            case 'spatialfilter'
                rowforspatialfilt = 2;
                rowforgrvote = 1;
            case 'gradientvote'
                rowforspatialfilt = 1;
                rowforgrvote = 2;
            case 'none'
                % see below -- still need neighborhood size, input 
                % as the first element in processparam, so color 'by hand'
                % below
                rowforspatialfilt = 1;
                rowforgrvote = 1;
                set(hlockobjsize, 'Value', true);  % force (irrelevant) 
                    % nsize to be same as bpfiltsize, which is used as the
                    % neighborhood size
                lockobjsize_Callback;
        end
        set(hbpfiltsizetext,'BackgroundColor', colorarray(rowforspatialfilt, :));
        set(hbpfiltsize,'BackgroundColor',colorarray(rowforspatialfilt, :));
        set(hnsizetext,'BackgroundColor',colorarray(rowforspatialfilt, :));
        set(hnsize,'BackgroundColor',colorarray(rowforspatialfilt, :));
        set(hgrobjsizetext,'BackgroundColor',colorarray(rowforgrvote, :));
        set(hgrobjsize,'BackgroundColor',colorarray(rowforgrvote, :));
        set(hgraddiropt,'BackgroundColor',colorarray(rowforgrvote, :));
        set(hgrthreshtext,'BackgroundColor',colorarray(rowforgrvote, :));
        if strcmp(processopt, 'none')
            set(hbpfiltsizetext,'BackgroundColor', colorarray(2,:));
            set(hbpfiltsize,'BackgroundColor',colorarray(2,:));
        end
        % set processing parameters (processparam variable)
        setprocessparam
        % update the processed images, if these are being displayed
        updateprocessandthresh;
    end

    function setprocessparam(source,eventdata)
        % set processing parameters (either for spatial filtering or
        % gradient voting), reading the appropriate values from the buttons
        % Note that the button callbacks are responsible for enforcing
        % consistency of the variables (e.g. slider and text box values
        % being the same)
        % Has input arguments, just to allow being used as a callback by
        % the hgraddiropt and hgrthreshtext buttons; inputs are not used.

        processopt = char(processoptarray(get(hprocessopt,'Value')));
        switch processopt
            case 'spatialfilter'
                % Spatial filtering
                processparam = [get(hbpfiltsize, 'Value') get(hnsize, 'Value')];
            case 'gradientvote'
                % Gradient voting
                graddiroptval = get(hgraddiropt, 'Value');
                switch graddiroptval
                    case 1
                        % both directions
                        graddir = 0;
                    case 2
                        % positive gradients
                        graddir = 1;
                    case 3
                        % negative gradients
                        graddir = -1;
                    otherwise
                        errordlg('Error -- bad graddiroptval')
                end
                processparam = [get(hgrobjsize, 'Value') graddir ...
                    str2double(get(hgrthreshtext, 'String'))];
            case 'none'
                % none, but set region size as first (only) element of
                % processparam
                processparam = get(hbpfiltsize, 'Value');
        end
    end

    function bpfiltsize_Callback(hObject, eventdata)
        % bandpass filter size -- called for both the text entry and the
        % slider.  Makes both button values equivalent, also updates 
        % neighborhood size if "locked," then calls
        % setprocessparam to update the processparam variable
        if hObject == hbpfiltsizetext
            % the text extry box has been called
            bpfiltsize = round(str2double(get(hObject,'string')));
            % Also update the slider:
            set(hbpfiltsize, 'Value', bpfiltsize);
        elseif hObject == hbpfiltsize
            % the slider has been called
            bpfiltsize = round(get(hObject,'value'));
            % Also update text entry box:
            set(hbpfiltsizetext, 'String', sprintf('%d', bpfiltsize));
            % to ensure an integer value, update the slider also!
            set(hbpfiltsize, 'Value', bpfiltsize);
        else
            errordlg('Error!  bad source for bpfiltsize_Callback!')
        end
        % if filter and neighborhood values are locked together, update
        % nsize buttons with bpfiltsize values
        if lockobjsize
            set(hnsize, 'Value', bpfiltsize);
            set(hnsizetext, 'String', sprintf('%d', bpfiltsize));
        end
        % set the processing parameter array with these values.
        setprocessparam
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end

    function nsize_Callback(hObject, eventdata)
        % neighborhood filter size -- called for both the text entry and the
        % slider.  Makes both button values equivalent, also updates 
        % filtering object size if "locked," then calls
        % setprocessparam to update the processparam variable
        if hObject == hnsizetext
            % the text extry box has been called
            nsize = round(str2double(get(hObject,'string')));
            % Also update the slider:
            set(hnsize, 'Value', nsize);
        elseif hObject == hnsize
            % the slider has been called
            nsize = round(get(hObject,'value'));
            % Also update text entry box:
            set(hnsizetext, 'String', sprintf('%d', nsize));
            % to ensure an integer value, update the slider also!
            set(hnsize, 'Value', nsize);
        else
            errordlg('Error!  bad source for nsize_Callback!')
        end
        % if filter and neighborhood values are locked together, update
        % nsize buttons with bpfiltsize values
        if lockobjsize
            set(hbpfiltsize, 'Value', nsize);
            set(hbpfiltsizetext, 'String', sprintf('%d', nsize));
        end
        % set the processing parameter array with these values.
        setprocessparam
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end

    function lockobjsize_Callback(hObject, source,eventdata)
        % make the filtering and neighborhood object sizes the same
        lockobjsize = get(hlockobjsize, 'Value');
        if lockobjsize
            set(hlockobjsize,'BackgroundColor', [0.8 1.0 0.6]);
            nsize = round(get(hbpfiltsize, 'Value'));
            set(hnsize, 'Value', nsize);
            set(hnsizetext, 'String', sprintf('%d', nsize));
            % update the filtered images, if these are being displayed
            updateprocessandthresh;
        else
            set(hlockobjsize,'BackgroundColor', [1.0 0.5 0.3]);
        end
    end

    function grobjsize_Callback(hObject, eventdata)
        % gradient voting object  size -- called for both the text entry and the
        % slider.  Makes both button values equivalent, then calls
        % setprocessparam to update the processparam variable
        if hObject == hgrobjsizetext
            % the text extry box has been called
            grobjsize = round(str2double(get(hObject,'string')));
            % Also update the slider:
            set(hgrobjsize, 'Value', grobjsize);
        elseif hObject == hgrobjsize
            % the slider has been called
            grobjsize = round(get(hObject,'value'));
            % Also update text entry box:
            set(hgrobjsizetext, 'String', sprintf('%d', grobjsize));
            % to ensure an integer value, update the slider also!
            set(hgrobjsize, 'Value', grobjsize);
        else
            errordlg('Error!  bad source for grobjsize_Callback!')
        end
        % set the processing parameter array with these values.
        setprocessparam
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end


    function whichthresh
        % determine which threshold value to use, based on thresholding
        % option
        threshopt = get(hthreshopt,'Value');
        switch threshopt
            case 1
                thresh = get(hthresh1,'Value');
            case 2
                thresh = -1.0*get(hthresh2,'Value');
                % fo5_rp.m interprets negative thresholds as "option 2"
                % inputs.
            case 3
                thresh = round(get(hthresh3,'Value'));
        end
    end
        
    function threshopt_Callback(source,eventdata)
        % set thresholding option
        threshopt = get(hthreshopt,'Value');
        ltgray = 0.95*[1 1 1];
        highlightcolor = [1 1 0.7];
        set(hthreshtext1,'BackgroundColor',ltgray);
        set(hthresh1,'BackgroundColor',ltgray);
        set(hthreshtext2,'BackgroundColor',ltgray);
        set(hthresh2,'BackgroundColor',ltgray);
        set(hthreshtext3,'BackgroundColor',ltgray);
        set(hthresh3,'BackgroundColor',ltgray);
        switch threshopt
            case 1
                set(hthreshtext1,'BackgroundColor',highlightcolor);
                set(hthresh1,'BackgroundColor',highlightcolor);
            case 2
                set(hthreshtext2,'BackgroundColor',highlightcolor);
                set(hthresh2,'BackgroundColor',highlightcolor);
            case 3
                set(hthreshtext3,'BackgroundColor',highlightcolor);
                set(hthresh3,'BackgroundColor',highlightcolor);
        end
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end

    function threshtext1_Callback(source,eventdata)
        % set thresh (text entry)
        thresh1 = str2double(get(source,'string'));
        set(hthresh1, 'Value', thresh1);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end
    function threshtext2_Callback(source,eventdata)
        % set thresh (text entry)
        thresh2 = str2double(get(source,'string'));
        set(hthresh2, 'Value', thresh2);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end
    function threshtext3_Callback(source,eventdata)
        % set thresh (text entry)
        thresh3 = round(str2double(get(source,'string')));
        set(hthresh3, 'Value', thresh3);
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end

    function thresh1_Callback(hObject, source,eventdata)
        thresh1 = get(hObject,'Value');
        % Also update text entry box:
        set(hthreshtext1, 'String', sprintf('%.4f', thresh1));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end
    function thresh2_Callback(hObject, source,eventdata)
        thresh2 = get(hObject,'Value');
        % Also update text entry box:
        set(hthreshtext2, 'String', sprintf('%.4f', thresh2));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end
    function thresh3_Callback(hObject, source,eventdata)
        thresh3 = round(get(hObject,'Value'));
        % Also update text entry box:
        set(hthreshtext3, 'String', sprintf('%.d', thresh3));
        % update the threshold value to use
        whichthresh;
        % update the filtered images, if these are being displayed
        updateprocessandthresh;
    end

    function ctrfitstr_Callback(hObject, source,eventdata)
        ctrfitstr = char(ctrfitstrarray(get(hObject,'Value')));
        fitstr = {ctrfitstr; orientstr; sizestr};
    end

    function orientstr_Callback(hObject, source,eventdata)
        orientstr = char(orientstrarray(get(hObject,'Value')));
        fitstr = {ctrfitstr; orientstr; sizestr};
    end

    function sizestr_Callback(hObject, source,eventdata)
        sizestr = char(sizestrarray(get(hObject,'Value')));
        fitstr = {ctrfitstr; orientstr; sizestr};
    end

    % Callbacks for Displaying processed and post-threshold images
    
    function processedA = calcprocessimg
        % calculate the 'processed' image, for neighborhood finding
        % either spatially filtered, or gradient voting option
        % Note that this simply copies the calculation in fo5_rp.m -- 
        %   recalculated when fo5_rp.m is called.
        switch processopt
            case 'spatialfilter'
                if processparam(1)>0
                    set(hmsg, 'String', 'Calling bpass.m to filter');
                    processedA = bpass(A,1,processparam(1));
                else
                    processedA = A;  % if filter parameter = 0, don't filter.
                end
            case 'gradientvote'
                set(hmsg, 'String', 'Calling gradientvote.m to process');
                processedA = gradientvote(A, processparam(1), processparam(2), processparam(3));
            case 'none'
                set(hmsg, 'String', 'No processing');
                processedA = A;
        end
    end
        
    function threshprocessA = calcthreshprocess
        % calculate the dilated image of local maxima that pass the
        % threshold.  
        % Calculate processed image even if previously calculated, in case
        % processing parameters have changed.
        processedA = calcprocessimg;  % uses bpfiltsize for filtering size
        % neighborhod size
        switch processopt
            case 'spatialfilter'
                nsize = processparam(2);
            case 'gradientvote'
                nsize = processparam(1);
            case 'none'
                nsize = processparam(2);
        end
        [y, x] = calcthreshpts(processedA, get(hthreshopt,'Value'), thresh, nsize);
        threshprocessAmask = false(size(processedA));
        threshprocessAmask(sub2ind(size(processedA), round(y), round(x)))=true;
        threshprocessAmask = imdilate(threshprocessAmask, strel('disk', floor(nsize/2)));
        threshprocessA = processedA.*threshprocessAmask;
    end

    function markedA = calcmarkers
        %calculate foreground and background markers
        %bwA = calcthreshprocess;
        %[ L, bgm, fgm4 ] = DomainWatershed( fixedA, 2, 3 );
        if islinkdone
            objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
        else
            objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
        end
        [ L, fgm, bgm ] = DomainWatershed( double(A), objsthisframe );
        markedA = A;
        markedA(L==0 | bgm | fgm) = 50;
        %markedA(imdilate(L==0,ones(3,3)) | fgm4) = 255;
        %markedA(imdilate(L==0,ones(3,3))) = 255;
        %markedA(fgm4) = 255;
        %markedA(bgm | fgm4 ) = 255;
        %         markedA(bgm) = 500;
    end

    function bilateralA = calcbilateralfilter
        if islinkdone
            objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
        else
            objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
        end
        A = img/max(max(A));
        B = bfilter2(A, 5, [3, .1]); %default values from bpfilter2
        bw = im2bw(B,graythresh(B));
        bilateralA(bw) = max(max(A));
    end

    function updateprocessandthresh
        % calls functions to recalculate filtered and post-threshold
        % images, and display
        if showprocess || showthreshprocess
            processedA = calcprocessimg;
        end
        if showthreshprocess
            threshprocessA = calcthreshprocess;
        end
%         if showwatershed
%             markedA = calcmarkers;
%         end
        if showbilateralfilter
            bilateralA = calcbilateralfilter;
        end
        displayimage
    end

    function dispprocess_Callback(hObject, source,eventdata)
        showprocess = get(hObject, 'Value');
        if showprocess
            processedA = calcprocessimg;  % calculate, even if processedA isn't empty, in case parameters have changed
            set(hdispthreshprocess, 'Value', false);  % turn off the other checkbox
%             set(hdispwater, 'Value', false);
        end
        displayimage;
    end

    function dispthreshprocess_Callback(hObject, source,eventdata)
        showthreshprocess = get(hObject, 'Value');
        % Could check to see if threshprocessA has already been calculated,
        % but should be careful that threshold parameter values haven't
        % changed.
        if showthreshprocess
            processedA = calcprocessimg;
            threshprocessA = calcthreshprocess;
            set(hdispprocess, 'Value', false);  % turn off the other checkbox
%             set(hdispwater, 'Value', false);
        end
        displayimage;
    end

    function dispwater_Callback(hObject, source, eventdata)
        showwatershed = get(hObject, 'Value');
        if showwatershed
            markedA = calcmarkers;
            set(hdispprocess, 'Value', false);
            set(hdispthreshprocess,'Value', false);
        end
        displayimage;
    end

    % Callback functions for Tracking (particle localization)
             
    function [] = trackthisfr_Callback(hObject, eventdata, handles)
        % Track a single frame (find objects)
        set(hmsg, 'String', 'Tracking started...'); 
        pause(0.1);  % does this help the display issue?
        tmpobj = fo5_rp(A, processopt, processparam, thresh, fitstr, ...
            get(h1pernhood, 'value'), [], lsqoptions);
        tmpobj(5,:) = currframe-frmin+1;
        objs = [objs tmpobj];
        trackthisdone(currframe-frmin+1) = true;  % note that tracking has been done
        savedobjs = objs;  % saved, in case culling and un-culling are done.
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(hmsg, 'String', 'Tracking completed!');
    end

    function [] = trackthisdone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the "tracking done" 
        % checkbox, reset its value to whatever the true array value is.
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
    end

    function [] = trackallfr_Callback(hObject, eventdata, handles)
        % Track all frames
        set(hmsg, 'String', 'Tracking of all images started...'); 
        pause(0.02);  % does this help the display issue?
        disp('Tracking of all images started...')

        progtitle = sprintf('TrackingGUI: Tracking using fo5_rp.m ...  '); 
        if (Nframes > 1)
            progbar = waitbar(0, progtitle);  % will display progress
        end
        objs = [];
        oldA = A;
        oldcurrframe = currframe;
        for j = frmin:frmax
            % Loop through all frames
            % similar code as in im2obj_rp.m
            if get(husenhoodctrs, 'value') && (j>frmin)
                % use x, y positions from
                % previous frame for neighborhood centers of this frame
                nhoodctrs = [tmpobj(1,:)' tmpobj(2,:)'];  
            else
                nhoodctrs = [];
            end
            if imagesloaded
                % all images already loaded
                tmpobj = fo5_rp(im(:,:,j), processopt, processparam, ...
                                thresh, fitstr, get(h1pernhood, 'value'), ...
                                nhoodctrs, lsqoptions);
            else
                % read from file; % use the variablea "A" and "currframe"
                currframe = j;
                loadimages;
                tmpobj = fo5_rp(A, processopt, processparam, ...
                                thresh, fitstr, get(h1pernhood, 'value'), ...
                                nhoodctrs, lsqoptions);
            end
            if ~isempty(tmpobj)
                tmpobj(5,:) = j-frmin+1;
            end
            objs = [objs tmpobj];
            % show progress -- not called if just one frame
            if mod(j-frmin+1,10)==0
                waitbar((j-frmin+1)/Nframes, progbar, ...
                    strcat(progtitle, sprintf('frame %d of %d', (j-frmin+1), Nframes)));
            end
            trackthisdone(j-frmin+1) = true;  % note that tracking has been done
        end
        A = oldA;
        currframe = oldcurrframe;
        if Nframes>1
            close(progbar)
        end
        savedobjs = objs;  % saved, in case culling and un-culling are done.
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
        set(hmsg, 'String', 'Tracking of all frames completed!');
        disp('Tracking of all frames completed!')
    end

    function [] = trackalldone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the  
        % checkbox, reset its value to what it should be
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
    end

    function [] = cleartracking_Callback(hObject, eventdata, handles)
        % Clear all the tracking output, and the linking output
        objs = [];
        objs_link = [];
        trackthisdone = false(Nframes,1);
        islinkdone = false;
        set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
        set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
        set(hlinkdone, 'Value', false);
        set(hmsg, 'String', 'Cleared tracking output for all frames');
    end

% Parameters for linking

    function  [] = link_Callback(hObject, eventdata, handles)
        
        % Linking objects into trajectories, using nnlink_rp
        set(hmsg, 'String', 'Linking tracks ...');
        disp('Linking tracks ...')
        
        if sum(trackthisdone)==length(trackthisdone)
            objs_link = nnlink_rp(objs, linkstep, linkmem, true);
            set(hmsg, 'String', 'Linking done ...');
            islinkdone = true;
            savedobjs_link = objs_link;  % saved, in case culling and un-culling are done.
        else
            % All frames need to be tracked before linking!
            set(hmsg, 'String', 'All frames need to be tracked before linking!');
        end
        set(hlinkdone, 'Value', islinkdone);
    end

    function [] = linkdone_Callback(hObject, eventdata, handles)
        % Don't do anything if the user clicks the "link done" 
        % checkbox, reset its value to whatever the true array value is.
        set(hlinkdone, 'Value', islinkdone);
    end

    function linksteptext_Callback(source,eventdata)
        % set max step size for linkage (text entry)
        linkstep = round(str2double(get(source,'string')));
    end

    function linkmemorytext_Callback(source,eventdata)
        % set memory size for linkage (text entry)
        linkmem = round(str2double(get(source,'string')));
    end

    function [] = clearlink_Callback(hObject, eventdata, handles)
        % Clear the linkages; keep objs matrix
        objs_link = [];
        islinkdone = false;
        set(hlinkdone, 'Value', islinkdone);
        set(hmsg, 'String', 'Cleared linkage of tracks');
    end

    function [] = culllink_Callback(hObject, eventdata, handles)
        % Cull objects from the objs or objs_link matrix based on histogram of object
        % width and histogram of gradient-line-distance, using cullobjs.m
        % If linking is done, apply to both objs and objs_link arrays; if
        % not apply to objs only
        prompt = {strcat('cullobjects.m:  Enter cutoffs for sigma, dmin as # std. above median,', ...
                   'separated by a space; leave empty for defaults;', ...
                   '-1 for input based on histograms'), ...
            'culling method: "rect" for rectangular cut in sigma/dmin space, "diag" for diagonal'};
        dlg_title = 'Culling parameters option'; num_lines= 1;
        def     = {num2str(-1), 'diag'};  % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        cutoffs = str2num(char(answer(1)));
        culloption = char(answer(2));
        if islinkdone
            [objs sigmamaxstd dmaxstd] = cullobjects(objs, cutoffs, culloption);
            objs_link = cullobjects(objs_link, [sigmamaxstd dmaxstd], culloption);
        else
            objs = cullobjects(objs, cutoffs, culloption);
        end
        set(hmsg, 'String', 'Objects culled');
    end

    function [] = unculllink_Callback(hObject, eventdata, handles)
        % Undo culling by reverting to saved matrices
        objs = savedobjs;
        objs_link = savedobjs_link;
        set(hmsg, 'String', 'Culling undone');
    end

%% Save output, or load previous calc. values
    function [] = saveoutput_Callback(hObject, eventdata, handles)
        % Save results in a MAT file
        % Dialog box for name
        prompt = {'Output File Name (*include* ".mat")'};
        dlg_title = 'Save output'; num_lines= 1;
        if isempty(outfile)
            def     = {'TrackGUIoutput.mat'};  % default value
        else
            def     = {outfile};  % default value
        end
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        outfile = char(answer(1));
        if ~isempty(outfile)
            set(hmsg, 'String', 'Starting to save variables...');
            save(outfile, 'Nframes', 'objs', 'objs_link', 'threshopt', ...
                'thresh', 'processopt', 'processparam', ...
                'trackthisdone', 'islinkdone');  % Save these variables
            set(hmsg, 'String', 'Done saving variables.');
        end
    end

    function [] = savetxt_Callback(hObject, eventdata, handles)
        % Save results in a simple text file; positions only
        % Dialog box for name
        if islinkdone
            prompt = {'Output File Name for simple Text Output (*include* extension)'};
            dlg_title = 'Save output'; num_lines= 1;
            def     = {'TrackGUIoutputTXT.txt'};  % default value
            answer  = inputdlg(prompt,dlg_title,num_lines,def);
            outfile = char(answer(1));
            if ~isempty(outfile)
                set(hmsg, 'String', 'Starting to save variables...');
                fo = fopen(outfile, 'w');
                utrk = unique(objs_link(6,:));
                for j=1:length(utrk)
                    xj = objs_link(1,objs_link(6,:)==utrk(j));
                    yj = objs_link(2,objs_link(6,:)==utrk(j));
                    for k=1: length(xj)
                        fprintf(fo, '%.3f\t', xj(k));
                    end
                    fprintf(fo, '\n');
                    for k=1: length(xj)
                        fprintf(fo, '%.3f\t', yj(k));
                    end
                    fprintf(fo, '\n');
                end
                fclose(fo);
                set(hmsg, 'String', 'Done saving variables.');
            end
        else
            beep
            set(hmsg, 'String', 'Linking must precede simple text output.');
        end
    end

    function [] = loadoutput_Callback(hObject, eventdata, handles)
        % Load obj data from previously saved MAT file
        % Dialog box for name
        prompt = {'Input File Name (*include* ".mat")'};
        dlg_title = 'Load output'; num_lines= 1;
        if isempty(outfile)
            def     = {'TrackGUIoutput.mat'};  % default value
        else
            def     = {outfile};  % default value
        end
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        infile = char(answer(1));
        if ~isempty(infile)
            S = load(infile);  % Load variables into a structure
            % -- necessary since this is a nested function
            % if Nframes is not saved, figure this out from range of row 5
            % in object matrix
            if ~isfield(S,'Nfames')
                S.Nframes = max(S.objs(5,:))-min(S.objs(5,:))+1;
            end
            if S.Nframes ~= Nframes
                set(hmsg, 'String', 'ERROR!  Number of frames is different');
                warndlg('Warning: Number of frames is inconsistent.');
                pause(2)
            else
                % reassign the variables!
                objs = S.objs;
                objs_link = S.objs_link;
                threshopt = S.threshopt;
                thresh = S.thresh;
                if isfield(S, 'processopt')
                    processopt = S.processopt;
                    processparam = S.processparam;
                else
                    % older version -- only spatial filtering as a
                    % processing option
                    processopt = 'spatialfilter';  % spatial filtering
                    if isfield(S, 'bpfiltsize')
                        processparam = [S.bpfiltsize S.nsize];
                    else
                        % even older version: one "objsize" variable
                        % instead of two for bpfiltsize and nsize
                        processparam = S.objsize*[1 1];
                    end
                end
                if isfield(S, 'islinkdone')
                    islinkdone = S.islinkdone;
                else
                    % field doesn't exist; guess from row 6
                    islinkdone = min(objs_link(6,:))>0;
                end
                if isfield(S, 'trackthisdone')
                    trackthisdone = S.trackthisdone;
                else
                    % field doesn't exist; guess from row 4 for each frame
                    utr_temp = unique(objs(5,:));
                    trackthisdone = false(Nframes,1);
                    for kk = 1:length(utr_temp)
                        trackthisdone(kk) = min(objs(4,(objs(5,:)==utr_temp(kk))))>0;
                    end
                end
                % update the filtered images, if these are being displayed
                updateprocessandthresh;
                set(hmsg, 'String', 'Done loading variables.');
                % set sliders and checkboxes
                set(hprocessopt, 'Value', find(strcmp(processoptarray,processopt)));
                switch processopt
                    case 'spatialfilter'
                        % Spatial filtering
                        set(hbpfiltsize, 'Value', processparam(1));
                        set(hbpfiltsizetext, 'String', sprintf('%d', processparam(1)));
                        set(hnsize, 'Value', processparam(2));
                        set(hnsizetext, 'String', sprintf('%d', processparam(2)));
                    case 'gradientvote'
                        % Gradient voting
                        set(hgrobjsize, 'Value', processparam(1));
                        set(hgrobjsizetext, 'String', sprintf('%d', processparam(1)));
                        set(hgraddiropt, 'Value', processparam(2));
                        set(hgrthreshtext, 'String', sprintf('%.3f', processparam(3)));
                    case 'none'
                        set(hbpfiltsize, 'Value', processparam(1));
                        set(hbpfiltsizetext, 'String', sprintf('%d', processparam(1)));
                end
                % Call processopt_Callback, just to set highlights color of
                % buttons
                processopt_Callback
                set(htrackthisdone, 'Value', trackthisdone(currframe-frmin+1));
                switch threshopt
                    case 1
                        set(hthresh1,'Value', thresh);
                        set(hthreshtext1, 'String', sprintf('%.4f', thresh));
                    case 2
                        set(hthresh2,'Value',-thresh);
                        set(hthreshtext2, 'String', sprintf('%.4f', -thresh));
                    case 3
                        set(hthresh3,'Value',thresh);
                        set(hthreshtext3, 'String', sprintf('%d', thresh));
                end
                set(hthreshopt, 'Value', threshopt);
                threshopt_Callback; % set the shading, etc.; a bit redundant
                set(htrackalldone, 'Value', sum(trackthisdone)==length(trackthisdone));
                set(hlinkdone, 'Value', islinkdone);  % update linkage variable
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call back functions for image plotting

    function displayimage(hObject, eventdata, handles)
        % Display the present frame
        if exist('im1', 'var')
            delete(im1); clear im1
        end
        if showthreshprocess
            im1 = imshow(threshprocessA, [], 'Parent', axeshandle);
        elseif showprocess
            im1 = imshow(processedA, [], 'Parent', axeshandle);
%         elseif showwatershed
%             im1 = imshow(markedA, [], 'Parent', axeshandle);
        elseif showbilateralfilter
            im1 = imshow(bilateralA, [], 'Parenst', axeshandle);
        else
            im1 = imshow(A, [], 'Parent', axeshandle);
        end
        set(axeshandle, 'Visible', 'off');
        colormap('gray');
        if ~isempty(displayRange)
            % Adjust display range
            caxis(displayRange*(2^bitDepth-1))
        end
        hold on
        
        % display track information on the image, if desired
        if trackthisdone(currframe-frmin+1)
            % tracking has been done
            if islinkdone
                objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
            else
                objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
            end
            if get(hdispcircles, 'Value')
                % plot circles on top, if tracking is done
                plot(objsthisframe(1,:), objsthisframe(2,:), 'o', 'color', [0.3 0.7 0.5])
                % If we're determining orientation, draw a line on each
                % object of length = major axis length and orientation =
                % found orientation
                if ~strcmp(orientstr, 'none')
                    % make an array of line starting and ending points
                    startx = objsthisframe(1,:) - objsthisframe(8,:).*cos(objsthisframe(7,:));
                    starty = objsthisframe(2,:) - objsthisframe(8,:).*sin(objsthisframe(7,:));
                    endx = objsthisframe(1,:) + objsthisframe(8,:).*cos(objsthisframe(7,:));
                    endy = objsthisframe(2,:) + objsthisframe(8,:).*sin(objsthisframe(7,:));
                    % draw lines
                    plot([startx; endx], [starty; endy], '-', 'color', [0.6 1.0 0.3])
                end
                if ~strcmp(sizestr, 'none')
                    %make array of lines of length equal to the object diameter
                    startx = objsthisframe(1,:) - objsthisframe(7,:);
                    starty = objsthisframe(2,:);
                    endx = objsthisframe(1,:) + objsthisframe(7,:);
                    endy = objsthisframe(2,:);
                    % draw lines
                    plot([startx; endx], [starty; endy], '-', 'color', [0.6 1.0 0.3])
                    for k = 1:length(objsthisframe(1,:))
                        if ~isnan(objsthisframe(7,k))
                            rectangle('Position',[objsthisframe(1,k)-objsthisframe(7,k) objsthisframe(2,k)-objsthisframe(7,k) 2*objsthisframe(7,k) 2*objsthisframe(7,k)], 'Curvature', [1,1],'EdgeColor', [0.6 1.0 0.3]);
                        end
                    end
                end
%                 (I think Tristan wrote this to superimpose a BW image of regions)
%                 if strcmp(sizestr, 'bilateral_filter')
%                     im1 = im1/max(max(im1));
%                     B = bfilter2(im1, 2*processparam(1)/2, [processparam(1)/2, 10]);
%                     bw = im2bw(B,graythresh(B));
%                     im1(bw) = max(max(im1));
% %                     markedA = A;
% %                     markedA(L==0 | bgm | fgm) = 50;
%                 end
            end
            if get(hdispIDs, 'Value')
                % show IDs
                if islinkdone
                    for j=1:size(objsthisframe,2)
                        text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(6,j)), 'color', [1.0 0.4 0.1])
                    end
                else
                    for j=1:size(objsthisframe,2)
                        text(objsthisframe(1,j), objsthisframe(2,j),num2str(objsthisframe(4,j)), 'color', [1.0 0.4 0.1])
                    end
                end
            end
            if get(hdisptracks, 'Value')
                % plot lines corresponding to the tracks of each object
                if islinkdone
                    % necessary
                    for j=1:size(objsthisframe,2)
                        allx = objs_link(1,objs_link(6,:)==objsthisframe(6,j));
                        ally = objs_link(2,objs_link(6,:)==objsthisframe(6,j));
                        plot(allx, ally, '-', 'color', [1 0.3 0.5])
                    end
                end
            end
        end
    end

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function output_txt = ObjectInfo_Callback(obj,eventdata)
        % Fill in the Object Properties table.
        % obj          Currently not used (empty)
        % event_obj    Object containing event data
        % output_txt   Data cursor text (string or cell array of strings)
        pos = get(eventdata,'Position');
        
        % Get the tracking information of the selected object
        if islinkdone
            objsthisframe = objs_link(:,objs_link(5,:)==currframe-frmin+1);
        else
            objsthisframe = objs(:,objs(5,:)==currframe-frmin+1);
        end
        % Euclidean distance matrix -- calculation method from Roland
        % Bunschoten's distance.m on the File Exchange
        aa=sum(pos'.*pos',1); 
        bb=sum([objsthisframe(1,:); objsthisframe(2,:)].*[objsthisframe(1,:); objsthisframe(2,:)],1);  
        d = sqrt(abs(aa( ones(size(bb,2),1), :)' + bb( ones(size(aa,2),1), :) - 2*pos*[objsthisframe(1,:); objsthisframe(2,:)]));
        % [avoid an extra function call] d = distance(pos', [objsthisframe(1,:); objsthisframe(2,:)]);
        [mind, imind] = min(d);
        thisobj = objsthisframe(:,imind);
        if size(thisobj,1)>1
            % not likely, but two closest objects; take the first one (arbitrary)
            thisobj = thisobj(:,1);
        end
        
        % Get the track or particle ID of the selected object
        
        % Properties of this object in this frame
        % When clicking on the object the only information in output_txt is 
        % shown on the figure
        % The rest of the information is sent to the data table.
        if islinkdone
            % Linkage across frames is done
            output_txt = thisobj(6);
        else
            % Not done, so use particle ID
            output_txt = thisobj(4);
        end
        set(t, 'Data', thisobj);  % for data table
    end



%Callback for the button group

    function selcbk(hObject,eventdata)
        switch get(eventdata.NewValue,'String') % Get Tag of selected object.
            case 'Zoom'
                % enable zoom feature.
                zoom(fGUI, 'on');
            case 'Object Information'
                %Enabling data cursor mode-This will allow us to collect data about
                %individual cells in the image.
                zoom(fGUI, 'off');
                set(dcm_obj,'Enable','on');
            case 'Pan'
                % enable pan feature
                pan(fGUI, 'on');
                % Continue with more cases as necessary. Might want to add a "grab'
                % feature to the image.
            otherwise
                % Code for when there is no match.
        end
        
    end

 end
