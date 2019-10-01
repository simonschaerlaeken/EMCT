function [ tsvFile ] = emcMidpoint( tsvFile, cfg)
% Computes the middle point between two markers of a list. Middle point is
% added to the marker list in the tsvFile
% 
% syntax
% tsvFile = emcMidpoint(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.transformMarker: cell array of markernames on which the baricenter and the midpoint will be computed
%     *.transformName: string containing the name of the group of markers creating the baricenter and middle point
%     *.betweenParameter: num defining the percetenge from the first point
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.transformMarker = {'a','b'};
% cfg.transformName = 'abc';
% cfg.betweenParameter = 80;
% tsvFile = emcMidpoint(tsvFile, cfg);
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
% Check param
if ~isfield(cfg, 'betweenParameter')
    cfg.betweenParameter = 50;
end
% addedMarker - keep in memory all the added markers - so that it can be
% checked that the tsvFile has been modified
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'addedMarker')
    tsvFile.info.addedMarker = {};
end
%% COMPUTATION AREA
%Compute Midle Point
markerIdxList = findIndexList(tsvFile.markerName, cfg.transformMarker);
tsvFileMarker = mcgetmarker(tsvFile,markerIdxList);
translation = tsvFileMarker.data(:,4:6)-tsvFileMarker.data(:,1:3);
middlepointData = tsvFileMarker.data(:,1:3)+translation*cfg.betweenParameter/100;
% Assign to middlepoint to TSV
tsvFile.data = horzcat(tsvFile.data, middlepointData);
tsvFile.nMarkers = tsvFile.nMarkers + 1;
markerName2 = cfg.transformName;
tsvFile.markerName{end+1} = markerName2; % Mid = Middle point
tsvFile.info.addedMarker{end+1} = markerName2;

end

