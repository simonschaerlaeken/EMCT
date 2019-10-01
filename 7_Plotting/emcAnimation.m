function [] = emcAnimation(tsvFile, cfg)
% Creates an animation of the motion capture data with the function
% mcanimate from the Motion Capture Toolbox, Copyright 2008, University of Jyvaskyla, Finland
% 
% syntax
% emcAnimation(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.connectionMatrix: list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     *.dataMapar: str path pointing at mcdemodata.mat from the Motion
%     Capture Toolbox, Copyright 2008, University of Jyvaskyla, Finland
%     [OPTIONAL]
%     *.fps: int number of frame per second at which the points will be displayed in the animation (default: 30)
%     *.msize: int size at which the points will be displayed (default: 2)
%     lines
%     
% output
% -
% 
% examples
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% cfg.dataMapar = 'C:/test/mcdemodata.mat';
% cfg.fps = 40;
% cfg.msize = 3;
% emcAnimation(tsvFile, cfg);
% 
% comments
% -
% 
% see also
% mcanimate
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

% connectionMatrix - connection matrix specifying all the connection
% between the markers
errorIfNotField(cfg, 'connectionMatrix');

% fps - fps at which the points will be displayed in the animation
if ~isfield(cfg, 'fps')||isempty(cfg.fps)
    disp('Since no cfg.fps was set up, arbirarily defined to 30.')
    cfg.fps = 30;
end
% msize - size at which the points will be displayed
if ~isfield(cfg, 'msize')||isempty(cfg.msize)
    disp('Since no cfg.msize was set up, arbirarily defined to 2.')
    cfg.msize = 3;
end
% dirMapar - path to mapar file in Toolbox Mocap - Necessary to run
if ~isfield(cfg, 'dataMapar')||~exist(cfg.dataMapar)
    warning('Please set up cfg.dirMapar to the directory containing the file "mapar" in the toolbox.')
    cfg.dataMapar = 'C:\Users\SCHAERLA\Google_Drive\PhD\6_Toolbox\61_Matlab\611_Motion_Capture_toolbox\MocapToolbox_v1.5\mocaptoolbox\mcdemodata.mat';
end


%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - Animation - fps: ' num2str(cfg.fps)]);
load(cfg.dataMapar, 'mapar') % Need to be in arg
mapar.msize = cfg.msize;
mapar.fps = cfg.fps;
mapar.conn = cfg.connectionMatrix;
mapar.scrsize = [1280,720];
mapar.cwidth = 2;
mapar.el = -15;
mapar.output = tsvFile.info.filename(1:end-4);
mcanimate(tsvFile, mapar);



end

