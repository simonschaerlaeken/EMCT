function [ tsvFile ] = emcFluidity(tsvFile, cfg)
% Computes the fluidity
% The principle of the curvature can be applied to the movement velocity 
% and its variation in time to calculate the movement fluidity (FI): high
% curvature of the speed trajectory in timemeans low fluidity, while 
% low curvaturemeans high fluidity. To calculate the fluidity we compute 
% the curvature of the tangential velocity of the desired joint. [Piana, 2013]
% 
% syntax
% tsvFile = emcFluidity(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.display = true
% cfg.displayUnit = 's'
% tsvFile = emcFluidity(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.fluidity
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end
% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
%% COMPUTATION AREA
disp('[PROCESSING] Fluidity')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
warning('off','all')
tsvFileSpeed = mctimeder(tsvFileMarker,1);
tsvFileAcc = mctimeder(tsvFileMarker,2);
warning('off','all')
fluidity = zeros(size(tsvFileMarker,1),tsvFileMarker.nMarkers);
for j = 1:tsvFileMarker.nMarkers
    tsvFileSpeedMarker = mcgetmarker(tsvFileSpeed, j);
    tsvFileAccMarker = mcgetmarker(tsvFileAcc, j);
    for i = 1:size(tsvFileMarker.data,1)
        fluidity(i,j) =  cross(tsvFileSpeedMarker.data(i,:),tsvFileAccMarker.data(i,:))/pow2(tsvFileSpeedMarker.data(i,:),3);
    end
end
% Average
fluidity = mean(fluidity,2);
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Fluidity'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, fluidity)
    else
        plot(fluidity)
    end
    title('Fluidity')
end

% Return the final value
tsvFile.processing.fluidity = fluidity;
end


