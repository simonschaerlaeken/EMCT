function [] = emcPlotBody(tsvFile, cfg)
%EMCPLOTBODY2D Function to plot the structure of the body in either 2D or
%3D
%   IN: tsvFile - structure defined as TSV
%       cfg     - configuration structure containing the options

%% CHECKING AREA

%% COMPUTATION AREA
% if no trace of previous projection to 2D
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'projected')||~tsvFile.info.projected
    emcPlotBody3D(tsvFile, cfg)
else
    emcPlotBody2D(tsvFile, cfg)
end

end

