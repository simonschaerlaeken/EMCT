function [ tsvFile ] = emcDirectness(tsvFile, cfg)
% Computes the directness - Movement Directness Index is computed from a
% trajectory drawn in the space by a joint as the ratio between the eu-
% clidean distance, calculated between the starting and the ending point of
% the trajectory, and the trajectory?s actual length. [Piana, 2013]
% 
% syntax
% tsvFile = emcDirectness(tsvFile, cfg)
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
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.display = true
% tsvFile = emcDirectness(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.directness
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



%% COMPUTATION AREA
disp('[PROCESSING] Directness')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
% Compute Directness
directness = 0;
for markerIdx = 1:length(tsvFileMarker.markerName)
    % Initialisation
    tsvFileMarkerTmp = mcgetmarker(tsvFileMarker,markerIdx);
    % If marker is full of NaN
    if isnan(tsvFileMarkerTmp.data)
        continue;
    end
    numerator = 0;
    denominator = 0;
    % Numerator computation according to formula
    for coordIdx = 1:3
        numerator = numerator + (tsvFileMarkerTmp.data(1,coordIdx)-tsvFileMarkerTmp.data(end,coordIdx)).^2;
    end
    numerator = sqrt(numerator);
    % Denominator computation according to formula
    for i=1:size(tsvFileMarkerTmp.data,1)-1
        tmp = 0;
        for coordIdx = 1:3
            tmp = tmp + (tsvFileMarkerTmp.data(i+1,coordIdx)-tsvFileMarkerTmp.data(i,coordIdx)).^2;
        end
        denominator = denominator + sqrt(tmp);
    end
    directness = directness + numerator/denominator;
end

% Return the final value
tsvFile.processing.directness = directness;
end


