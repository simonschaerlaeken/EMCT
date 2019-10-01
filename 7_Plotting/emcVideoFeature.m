function [ ] = emcVideoFeature( tsvFile, cfg )
% Plots an animation of both the body moving and a defined feature evolving over
% time
% 
% syntax
% emcVideoFeature(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.connectionMatrix: list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     *.videoFeature: str feature name to be displayed
%     [OPTIONAL]
%     *.step: int number of frames to skip in between two timestamps (default: 50)
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
% emcVideoFeature(tsvFile, cfg);
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
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end

% videoFeature - define the feature to use
errorIfNotField(cfg, 'videoFeature')
errorIfNotField(tsvFile, 'processing')
errorIfNotField(tsvFile.processing, cfg.videoFeature)

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
if ~isfield(tsvFile, 'videoFrameRate')
    disp('videoFrameRate not defined - use default = 15')
    cfg.videoFrameRate = 15; % 'body feature'
end

%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - VIDEO']);

F(length(1:cfg.step:size(tsvFile.data,1))) = struct('cdata',[],'colormap',[]);
title = tsvFile.info.filename;
title = [title, ' - ', cfg.videoFeature];
% Find max and min for axis limits for mocap data;
% X
xMin = min(min(tsvFile.data(:,1:3:size(tsvFile.data, 2))));
xMax = max(max(tsvFile.data(:,1:3:size(tsvFile.data, 2))));
% Y
yMin = min(min(tsvFile.data(:,2:3:size(tsvFile.data, 2))));
yMax = max(max(tsvFile.data(:,2:3:size(tsvFile.data, 2))));
% Z
zMin = min(min(tsvFile.data(:,3:3:size(tsvFile.data, 2))));
zMax = max(max(tsvFile.data(:,3:3:size(tsvFile.data, 2))));
% Save it in cfg
cfg.limitAxis = [xMin, xMax, yMin, yMax, zMin, zMax];
    
figure('Name',title,'NumberTitle','off','units','normalized','outerposition',[0 0 0.5 1])
h1 = subplot(2,1,1);
h2 = subplot(2,1,2);
idx = 1;
for frame = 1:cfg.step:size(tsvFile.data,1),
    axesPlotBodyInstant(h1, tsvFile, cfg, frame);
    axesPlotFeatureInstant(h2, tsvFile.processing.(cfg.videoFeature), frame);
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

