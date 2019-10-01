function [ tsvFile ] = emcAlignMarker( tsvFile, cfg )
% Align a marker dimension onto others markers
% 
% syntax
% emcPlotBody3D(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.markerToAlign: string containing the name of the marker to align
%     *.markerToAlignTo: cell of strings containing the marker(s) to align
%     *.dimensionToAlign: num defining the dimension to align (1, 2, or 3 corresponding to X, Y, or Z)
%     it to
%     [OPTIONAL]
%     -
%     
% output
% -
% 
% examples
% cfg.markerToAlign = 'a'
% cfg.markerToAlignTo = {'b','c'};
% cfg.dimensionToAlign = 3;
% emcAlignMarker(tsvFile, cfg);
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
errorIfNotField(cfg, 'markerToAlign')
errorIfNotField(cfg, 'markerToAlignTo')
errorIfNotField(cfg, 'dimensionToAlign')
%% COMPUTATION AREA
% Isolate markers to align it to
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.markerToAlignTo);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
dataToAlignTo = centroid3D(tsvFileMarker.data);
dataToAlignTo = dataToAlignTo(:, cfg.dimensionToAlign); % Take only the dimension of interest
% Find the marker to align
featMarkerIndex = findIndexList(tsvFile.markerName, {cfg.markerToAlign});
indexToChange = (featMarkerIndex(1)-1)*3+cfg.dimensionToAlign;
tsvFile.data(:,indexToChange) = dataToAlignTo;
end

