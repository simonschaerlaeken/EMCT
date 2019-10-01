function [ ] = emcVideoFeatureGroup( varargin )
% Plots an animation of both the bodies of a group of people moving and a defined feature evolving over
% time
% 
% syntax
% emcVideoFeatureGroup(tsvFile1, tsvFile2, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.connectionMatrix: list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     *.videoFeature: str feature name to be displayed
%     [OPTIONAL]
%     *.fps: int number of frames to skip in between two timestamps (default: 50)
%     *.videoName: str name of the video to be saved
%     *.videoFrameRate: int framerate at which the video is recorded
%     
% output
% -
% 
% examples
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% cfg.videoFeature = 'speed';
% cfg.step = 30;
% cfg.videoName = 'testVideo.avi';
% cfg.videoFrameRate = 45;
% emcVideoFeature(tsvFile1, tsvFile2, cfg);
% 
% comments
% Feature must have been computed beforehand using funciton in Processing
% Directory
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if numel(varargin)<2,
    error('ERROR: No file or no cfg as input')
elseif numel(varargin)<3,
    disp('[INFO] Only one file to visualise, last argument is cfg')
end
cfg = varargin{end};
tsvFileIdxList = 1:numel(varargin)-1;
for tsvFileIdx = tsvFileIdxList,
    % filename - for title
    if ~isfield(varargin{tsvFileIdx}, 'info')||~isfield(varargin{tsvFileIdx}.info, 'filename'),
        filename = strsplit(varargin{tsvFileIdx}.filename, filesep);
        filename = filename{end};
        varargin{tsvFileIdx}.info.filename = filename;
    end
end

% videoFeature - define the feature to use
errorIfNotField(cfg, 'videoFeature')
for tsvFileIdx = tsvFileIdxList,
    errorIfNotField(varargin{tsvFileIdx}, 'processing')
    errorIfNotField(varargin{tsvFileIdx}.processing, cfg.videoFeature)
end

% connectionMatrix - table containing all the connexions (sticks) between
% markers. Necessary to run
errorIfNotField(cfg, 'connectionMatrix')

% step - define the number of frames to skip in between two timestamps.
if ~isfield(cfg, 'step')
    cfg.step = 50;
end
% videoName - name of the saved video.
if ~isfield(cfg, 'videoName')
    cfg.videoName = ['video_', cfg.videoFeature,'.avi'];
end
% videoFrameRate - framerate at which the video is recorded.
if ~isfield(cfg, 'videoFrameRate'),
    disp('videoFrameRate not defined - use default = 15')
    cfg.videoFrameRate = 15; % 'body feature'
end

%% COMPUTATION AREA
title = '';
for tsvFileIdx = tsvFileIdxList,
    title = [varargin{tsvFileIdx}.info.filename, ' '];
end
disp(['[PLOT] ', title ' - VIDEO']);

% TotalFrame taken on the first body
totalFrame = size(varargin{1}.data,1);

F(length(1:cfg.step:totalFrame)) = struct('cdata',[],'colormap',[]);
title = [title, ' - ', cfg.videoType];
title = [title, ' - ', cfg.videoFeature];

% Find max and min for axis limits for mocap data;
xMin = zeros(1,numel(tsvFileIdxList)); xMax = zeros(1,numel(tsvFileIdxList));
yMin = zeros(1,numel(tsvFileIdxList)); yMax = zeros(1,numel(tsvFileIdxList));
zMin = zeros(1,numel(tsvFileIdxList)); zMax = zeros(1,numel(tsvFileIdxList));
for tsvFileIdx = tsvFileIdxList,
    tsvFile = varargin{tsvFileIdx};
    % X
    xMin(tsvFileIdx) = min(min(tsvFile.data(:,1:3:size(tsvFile.data, 2))));
    xMax(tsvFileIdx) = max(max(tsvFile.data(:,1:3:size(tsvFile.data, 2))));
    % Y
    yMin(tsvFileIdx) = min(min(tsvFile.data(:,2:3:size(tsvFile.data, 2))));
    yMax(tsvFileIdx) = max(max(tsvFile.data(:,2:3:size(tsvFile.data, 2))));
    % Z
    zMin(tsvFileIdx) = min(min(tsvFile.data(:,3:3:size(tsvFile.data, 2))));
    zMax(tsvFileIdx) = max(max(tsvFile.data(:,3:3:size(tsvFile.data, 2))));
end
% Save it in cfg
cfg.limitAxis = [min(xMin), max(xMax), min(yMin), max(yMax), min(zMin), max(zMax)];

figure('Name',title,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
h1 = subplot(2,1,1);
h2 = subplot(2,1,2);
idx = 1;
for frame = 1:cfg.step:totalFrame,
    % Create all the tsv from VARARGIN
    inputTsvFile = '';
    for tsvFileIdx = tsvFileIdxList,
        eval(['tsvFile',num2str(tsvFileIdx),' = varargin{', num2str(tsvFileIdx),'};'])
        inputTsvFile = [inputTsvFile,'tsvFile',num2str(tsvFileIdx),','];
    end
    inputTsvFile = [inputTsvFile,' cfg, frame'];
    % Plot body
    eval(['axesPlotBodyInstant(h1, ', inputTsvFile, ');']);
    % Plot feature on the first Body since all of them have the same group
    % features
    axesPlotFeatureInstant(h2, varargin{1}.processing.(cfg.videoFeature), frame);
    F(idx) = getframe(gcf);
    idx = idx + 1;
end

myVideo = VideoWriter(cfg.videoName);
myVideo.FrameRate = cfg.videoFrameRate;  % Default 15
myVideo.Quality = 100; 
open(myVideo);
writeVideo(myVideo, F);
close(myVideo);

end

