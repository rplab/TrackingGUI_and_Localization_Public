% getnumfilelist.m
%
% Function to get the "base" file name and starting and ending frames for a
% series of images with filenames 'NAMExxxx.tif', where NAME is the same for
% all images and xxxx is a sequential numbering.
%
% Note: the filenames must end with the number (other than the extention)
%
% If the second file name is the same as the first, return empty arrays for
% the starting and ending frames. 
% Can also be used for multi-page TIFFs, in which case the second file name
% is output as empty.  Function detects that the "first" filname chosen is
% a multipage TIFF
% Allows frmax to be less than frmin -- whatever calls this should deal
% with this issue. However, determine "formatstr" based on the proper
% ordering.
%
% Originally extracted from TIFFseries.m
%
% Inputs:
%    FileName1, FileName2 : (Optional) The names of the starting and ending
%         image files. Include the extension. If input, bypass dialog box.
%         Leave FileName2 empty if importing a multi-page TIFF
%        
% Outputs
%    baseFileName : Base file name string  -- i.e. the part without numbers. Whole file name for a multipage TIFF
%    frmin : min frame number (i.e. min xxxx) -- "1" for a multipage TIFF
%    frmax : max frame number (i.e. min xxxx) -- "#frames" for a multipage TIFF
%    formatstr : formatting string to create file names
%    FileName1 : string of first file name 
%    FileName2 : string of last file name 
%    PathName1 : path of the (first) file
%    ext : file extension (e.g. '.tif')
%    ismultipage : boolean true if considering a multipage TIFF (redundant
%            with FileName2 being empty)
%
% EXAMPLE
% To use the output variables to create file names , e.g.:
%        k = 34;
%        framestr = sprintf(formatstr, k);
%        FileName = strcat(fbase, framestr, ext);
%        A  = imread(FileName);  % image
%
% EXAMPLE
% To use the output variables to load an image from a mulitpage TIFF:
%       k = 34;
%       A  = imread(strcat(fbase, ext), k);
% Raghuveer Parthasarathy
% Sept. 15, 2010
% July 10, 2011: Get rid of initial dialog, and go immediately to file 
%    selection box.Allow multi-page TIFF input
% May 12, 2020: Allow different file lengths (e.g. ABC_1.tif, ...,
%    ABC_12.tif)
% last modified: June 23, 2020

function [baseFileName, frmin, frmax, formatstr, FileName1, FileName2, PathName1, ext, ismultipage] = ...
    getnumfilelist(FileName1, FileName2)

firstdir = pwd;
if nargin > 0
    % Filename(s) are specified, avoid dialog box.
    [PathName1] = fileparts(FileName1);  % to get the extension (includes ".")
else
    % Dialog boxe for filenames
    [FileName1, PathName1] = uigetfile('*.*', 'First .tif image, or multi-page TIFF...');
    NumberOfFrames1 = numel(imfinfo(strcat(PathName1,FileName1)));
    if  NumberOfFrames1 > 1
        % Multipage TIFF
        FileName2 = [];
    else
        cd(PathName1)
        [FileName2] = uigetfile('*.*', 'Last .tif image of sequence...');
        cd(firstdir)
    end
end

% If Filename2 is empty, we must be importing a multi-page TIFF.
if isempty(FileName2)
    ismultipage = true;
else
    ismultipage = false;
end

%% Get file info
[~, FileName1String, ext] = fileparts(FileName1);  % to get the extension (includes ".")
if ismultipage
    frmin = 1;
    frmax = NumberOfFrames1;
    formatstr = [];
    baseFileName = FileName1String;
else
    % series of 2D TIFF images
    [~, ~, ext] = fileparts(FileName1);  % to get the extension (includes ".")
    if strcmp(FileName1, FileName2)
        % Both file names are the same -- user has selected only one file
        disp('Only one file selected! Leaving frmin, frmax empty.');
        frmin = [];
        frmax = [];
        formatstr = [];
        baseFileName = FileName1String;
    else
        % the usual case(starting and ending images)
        % Find the character position before the terminating numbers
        lastChar = findLastLetterPosition(FileName1String);
        if isempty(lastChar)
            lastChar = 0; % necessary if the filename consists only of numbers
        end
        frmin = str2num(FileName1String(lastChar+1: end)); %#ok<ST2NM>
        baseFileName = FileName1String(1:lastChar);
        % Ending frame
        [~, FileName2String] = fileparts(FileName2);  % to get the extension (includes ".")
        lastChar2 = findLastLetterPosition(FileName2String);
        if isempty(lastChar2)
            lastChar2 = 0; % necessary if the filename consists only of numbers
        end
        frmax = str2num(FileName2String(lastChar2+1: end)); %#ok<ST2NM>
        if ~strcmpi(baseFileName, FileName2String(1:lastChar2))
            % "base" file names don't match!
            errordlg('Error: "base" file names don''t match!')
            baseFileName = [];
            formatstr = [];
        else
            if (length(FileName1String)-lastChar) == (length(FileName2String)-lastChar2)
                % same number of digits in both files' numeric parts
                digits = length(FileName1String)-lastChar;
                formatstr = strcat('%0', num2str(digits), 'd');
            else
                formatstr = '%d';
            end
        end
    end
end

end

function lastChar = findLastLetterPosition(string)
    % Find the character position before the terminating numbers
    nChar = length(string);  % number of characters in file name
    numPositions = regexp(string,'\d');
    ismember(1:nChar, numPositions);
    lastChar = find(diff(ismember(1:nChar, numPositions))==1, 1, 'last');
end
