function [ tsvFile ] = emcAlignHeightMarkerBaricenter( tsvFile, cfg )
% Aligns markers to the height of the baricenter of a collection of markers.
% 
% syntax
% tsvFile = emcAlignHeightMarkerBaricenter(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.transformMarker: string containing the path to the output directory 
%     *.alignMarker: cell array of markers names to be deleted
%     [OPTIONAL]
%     -
%   
% output
% []
% 
% examples
% cfg.transformMarker = {'a','b'};
% cfg.alignMarker = {'c'};
% tsvFile = emcAlignHeightMarkerBaricenter(tsvFile, cfg);
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
% Check error
errorIfNotField(cfg, 'transformMarker') % Marker to compute the baricenter and the midpoint
if isstring(cfg.transformMarker) % Make it into a cell to be used
    cfg.transformMarker = {cfg.transformMarker};
end
% Check error
errorIfNotField(cfg, 'alignMarker') % Marker to compute the baricenter and the midpoint
if isstring(cfg.alignMarker) % Make it into a cell to be used
    cfg.alignMarker = {cfg.alignMarker};
end

%% COMPUTATION AREA
for markerIdx = 1:numel(cfg.transformMarker)
    % Find transform ind
    transformMarkerIdxList = findIndexList(tsvFile.markerName, cfg.transformMarker{markerIdx});
    % Compute baricenter
    alignMarkerIdx = findIndexList(tsvFile.markerName, cfg.alignMarker);
    tsvFileAlignMarker = mcgetmarker(tsvFile,transformMarkerIdxList);
    baricenter = centroid3D(tsvFileAlignMarker.data);
    % Replace the Z coordinate
    tsvFile.data(:,alignMarkerIdx*3) = baricenter(:,3);
end
end

