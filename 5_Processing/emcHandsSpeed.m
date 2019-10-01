function [ tsvFile ] = emcHandsSpeed(tsvFile, cfg)
% Computes the speed of the hands
% 
% syntax
% tsvFile = emcHandsSpeed(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.handsMarkers: cell array containing the markernames for the hands
%     *.euclidianFlag: boolean to determine if the script uses euclidian
%     norm after derivative
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.handsMarkers = {'a','b'};
% cfg.display = true
% cfg.displayUnit = 's'
% cfg.euclidianFlag = true;
% tsvFile = emcHandsSpeed(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.speed
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
% If euclidian norm or raw data as output
errorIfNotField(cfg, 'euclidianFlag') 

% Check if missing
errorIfNotField(cfg, 'handsMarkers')



%% COMPUTATION AREA
disp('[PROCESSING] Hands Speed')
% Set up cfg
cfg.deriv = 1;
cfg.featMarker = cfg.handsMarkers;
% Call derivative
tsvFile = emcDerivative(tsvFile, cfg);

end


