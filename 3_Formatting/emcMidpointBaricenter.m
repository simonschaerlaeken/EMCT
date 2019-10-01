function [ tsvFile ] = emcMidpointBaricenter( tsvFile, cfg)
% Computes the baricenter of a list of multiple markers and the middle point between
% the first two markers of the list. The baricenter and middle point are
% added to the marker list in the tsvFile
% 
% syntax
% tsvFile = emcMidpointBaricenter(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.transformMarker: cell array of markernames on which the baricenter and the midpoint will be computed
%     *.transformName: string containing the name of the group of markers creating the baricenter and middle point
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.transformMarker = {'a','b','c'};
% cfg.transformName = 'abc';
% tsvFile = emcMidpointBaricenter(tsvFile, cfg);
% 
% comments
% The markers created are a combination of the transformName and either "Bar" for Baricenter or "Mid" for Middle point
% 
% see also
% emcVerticalLine
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% Check error
% transformMarker: cell array of markernames on which the baricenter and the midpoint will be computed
errorIfNotField(cfg, 'transformMarker') 
% Check error
errorIfNotField(cfg, 'transformName') % Marker to compute the baricenter and the midpoint
% addedMarker - keep in memory all the added markers - so that it can be
% checked that the tsvFile has been modified
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'addedMarker')
    tsvFile.info.addedMarker = {};
end
%% COMPUTATION AREA

%Compute Baricentre
markerIdxList = findIndexList(tsvFile.markerName, cfg.transformMarker);
tsvFileMarker = mcgetmarker(tsvFile,markerIdxList);
baricentreData = centroid3D(tsvFileMarker.data);

%Compute Midle Point
markerIdxList = findIndexList(tsvFile.markerName, cfg.transformMarker(1:2));
tsvFileMarker = mcgetmarker(tsvFile,markerIdxList);
middlepointData = centroid3D(tsvFileMarker.data);

% Assign to baricentre to TSV
tsvFile.data = horzcat(tsvFile.data, baricentreData);
tsvFile.nMarkers = tsvFile.nMarkers + 1;
markerName1 = [cfg.transformName 'Bar'];
tsvFile.markerName{end+1} = markerName1; % Bar = Baricentre
tsvFile.info.addedMarker{end+1} = markerName1;

% Assign to middlepoint to TSV
tsvFile.data = horzcat(tsvFile.data, middlepointData);
tsvFile.nMarkers = tsvFile.nMarkers + 1;
markerName2 = [cfg.transformName 'Mid'];
tsvFile.markerName{end+1} = markerName2; % Mid = Middle point
tsvFile.info.addedMarker{end+1} = markerName2;


end

