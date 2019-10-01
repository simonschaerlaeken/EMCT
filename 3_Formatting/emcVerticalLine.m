function [ tsvFile ] = emcVerticalLine( tsvFile, cfg )
% Creates a second marker right above a specific marker in order to draw a vertical line between them - 
% 
% syntax
% tsvFile = emcVerticalLine( tsvFile, cfg );
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.transformMarker: cell array or string containing markername(s)
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.transformMarker = {'a'};
% tsvFile = emcMidpointBaricenter(tsvFile, cfg);
% 
% comments
% The markers created are a combination of the markername and "Vert" for Vertical line
% This can be further use to compare angle to a straight line or to have an idea
% of how bended the movement is compared to a straight vertical line
% 
% see also
% emcMidpointBaricenter
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% Check error
errorIfNotField(cfg, 'transformMarker') % Marker to compute the baricenter and the midpoint
if isstring(cfg.transformMarker) % Make it into a cell to be used
    cfg.transformMarker = {cfg.transformMarker};
end
% addedMarker - keep in memory all the added markers - so that it can be
% checked that the tsvFile has been modified
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'addedMarker')
    tsvFile.info.addedMarker = {};
end
%% COMPUTATION AREA
for markerIdx = 1:numel(cfg.transformMarker)
    % Get the x ,y and z coordinates of the marker
    markerIdxList = findIndexList(tsvFile.markerName, cfg.transformMarker{markerIdx});
    tsvFileMarker = mcgetmarker(tsvFile,markerIdxList);
    dataMarker = tsvFileMarker.data(:,:);
    % Define value for second point to create the line
    dataMarker(:,3) = dataMarker(:,3) + 1; % +1 to differentiate and stay on the Z axis
    tsvFile.data = horzcat(tsvFile.data, dataMarker);
    tsvFile.nMarkers = tsvFile.nMarkers + 1;
    markerName = [cfg.transformMarker{markerIdx}, 'Vert']; % Vert = vertical 
    tsvFile.markerName{end+1} = markerName;
    tsvFile.info.addedMarker{end+1} = markerName;
end

