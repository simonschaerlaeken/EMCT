function [] = emcPlotMarker(tsvFile, cfg)
% Plots the motion capture data of different markers selected
% in cfg in different colours on the three axis X, Y, and Z
% 
% syntax
% emcPlotMarker(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.plotMarkerList: cell array with markers to be plotted (if not defined all markers will be used)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     
% output
% -
% 
% examples
% cfg.plotMarkerList = {'a','b','c'};
% cfg.displayUnit = 's'
% emcPlotMarker(tsvFile, cfg);
% 
% comments
% -
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
% plotMarkerList - choice of marker to display
if ~isfield(cfg, 'plotMarkerList')||isempty(cfg.plotMarkerList)
    warning('No marker was define, by default all markers are printed. Define markers to be printed using "cfg.plotMarkerList"')
    cfg.plotMarkerList = tsvFile.markerName;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - Compare Markers']);
% Extract data using just the markers we want to compare
listMarkerSelectedIdx = findIndexList(tsvFile.markerName,cfg.plotMarkerList);
tsvFileMarker = mcgetmarker(tsvFile,listMarkerSelectedIdx);
% Plot every markers on a graph for the different parameters
hFig = figure('Name',[tsvFile.info.filename ' -- Compare Markers'], 'NumberTitle','off'); clf; hold on;
set(hFig, 'Position' , [230 230 1300 700])
% Create a different color for each value
colors_individual = distinguishable_colors(length(tsvFileMarker.markerName));

for i = 1:3:length(tsvFileMarker.markerName)*3
    % Print Graph
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        hold on,ax1 = subplot(3,1,1);plot(time, tsvFileMarker.data(:,i),'color' ,colors_individual((i-1)/3+1,:)),title('x')
        hold on,ax2 = subplot(3,1,2);plot(time, tsvFileMarker.data(:,i+1),'color' ,colors_individual((i-1)/3+1,:)),title('y')
        hold on,ax3 = subplot(3,1,3);plot(time, tsvFileMarker.data(:,i+2),'color' ,colors_individual((i-1)/3+1,:)),title('z')
    else
        hold on,ax1 = subplot(3,1,1);plot(tsvFileMarker.data(:,i),'color' ,colors_individual((i-1)/3+1,:)),title('x')
        hold on,ax2 = subplot(3,1,2);plot(tsvFileMarker.data(:,i+1),'color' ,colors_individual((i-1)/3+1,:)),title('y')
        hold on,ax3 = subplot(3,1,3);plot(tsvFileMarker.data(:,i+2),'color' ,colors_individual((i-1)/3+1,:)),title('z')
    end
    
end
legend(cfg.plotMarkerList)

% Set global Title
% annotation(hFig, 'textbox',[0.49 0.97 0.12 0.032],...
%     'String',title_fig,...
%     'FitBoxToText','on' ,...
%     'BackgroundColor',[1 1 1]);
end



