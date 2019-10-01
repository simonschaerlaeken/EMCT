function [ tsvFile ] = emcUpwardMovement(tsvFile, cfg)
% Computes the amount of upward movement, everytime a marker is going up on
% the Z axis, the value for the feature increases
% 
% syntax
% tsvFile = emcUpwardMovement(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% tsvFile = emcUpwardMovement(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.upwardMovement
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

%% COMPUTATION AREA
disp('[PROCESSING] Upward Movement')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
% focus on marker bow up (z ie vertical component)
%initiate
upwardMovement = zeros(3,length(cfg.featMarker));
% Obtain velocity
warning('off','all')
tsvFileSpeed = mctimeder(tsvFileMarker,1);
warning('on','all')
% Use vertical axis
% param_computation = param_computation*3;
for idxMarkerZ = 1:length(cfg.featMarker)
    markerZ = featMarkerIndex(idxMarkerZ)*3;
    tsvFileDataMarkerZ = tsvFileSpeed.data(:,markerZ);
    for i=1:length(tsvFileDataMarkerZ)
        if tsvFileDataMarkerZ(i)>0 % Positive speed in Z coordinate
           upwardMovement(1,idxMarkerZ)=upwardMovement(1,idxMarkerZ)+1;
        elseif tsvFileDataMarkerZ(i) < 0 % Negative speed in Z coordinate
           upwardMovement(2,idxMarkerZ)=upwardMovement(1,idxMarkerZ)+1; 
        else % No speed in Z coordinate 
           upwardMovement(3,idxMarkerZ)=upwardMovement(1,idxMarkerZ)+1;
        end
    end
end

% Return the final value
tsvFile.processing.upwardMovement = upwardMovement;
end


