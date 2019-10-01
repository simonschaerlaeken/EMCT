function [] = emcHeatmap( tsvFile, cfg )
% Computes the heatmap of the baricenter motion seen only on a 2D plane (x-y) from above. 
% 
% syntax
% tsvFile = emcHeatmap(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.heatmapFeature: str feature name to modulate the value of the heatmap
%     
% output
% -
% 
% examples
% cfg.heatmapFeature = 'speed';
% cfg.reduceFlag = true;
% cfg.reduceResolution = [2,3];
% emcHeatmap(tsvFile, cfg);
% 
% comments
% -
% 
% see also
% emcHeatmapData
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% Checking Area
if ~isfield(cfg, 'heatmapFeature')||isempty(cfg.heatmapFeature)
    cfg.heatmapFeature = '';
else
    if ~isfield(tsvFile.processing, cfg.heatmapFeature)
        error(['Feature ', cfg.heatmapFeature, 'not present for this tsvFile'])
    end
end
%% COMPUTATION AREA
if isempty(cfg.heatmapFeature)
    space = emcHeatmapData(tsvFile.data, cfg);
else
    space = emcHeatmapData(tsvFile.data, tsvFile.processing.(cfg.heatmapFeature), cfg);
end

