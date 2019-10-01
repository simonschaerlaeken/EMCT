function [ tsvFile ] = emcAbsolutePosition(tsvFile, cfg)
% Computes absolute position of tsvFile.
% Useful if need to be modify with for example preprocessing steps such
% as mccenter.
% 
% syntax
% tsvFile = emcAbsolutePosition(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.headMarkers: cell array with markers of the head
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.headMarkers = {'head1','head2'};
% tsvFile = emcAbsolutePosition(tsvFile, cfg);
% 
% comments
% Stores the absolute position in tsvFile.info.absolutePosition
% 
% see also
% emcLoadSingle
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% headMarkers - Defined the markers of the head. If exist and not empty,
% then compute a baricenter for head
if isfield(cfg, 'headMarkers') % Check argument for prefix
    if isempty(cfg.headMarkers)
        disp('Head Markers for absolute position are not defined')
    end
end
%% COMPUTATION AREA
% First data
if size(tsvFile.data,2) >= 3
    tsvFile.info.absolutePosition.firstMarker = tsvFile.data(:,1:3);
else
    disp('Data size is inferior than 3')
    tsvFile.info.absolutePosition.firstMarker = tsvFile.data;
end

% Baricenter
tsvFile.info.absolutePosition.baricenter = centroid3D(tsvFile.data);

% Head
if isfield(cfg, 'headMarkers')
    tsvFileHead = mcgetmarker(tsvFile, findIndexList(tsvFile.markerName, cfg.headMarkers));
    tsvFile.info.absolutePosition.head = centroid3D(tsvFileHead.data);
end
end

