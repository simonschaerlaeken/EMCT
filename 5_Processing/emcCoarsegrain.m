function [tsvFile] = emcCoarsegrain(tsvFile,cfg) 
% Loads a single tsv file from the filename into a structure
% 
% syntax
% tsvFile = emcCoarsegrain(tsvFile,cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.coarsegrainType: str type of computation of the data -- 'mean' /
%     'sum'. The computation is applied after the segmentation.
%     *.scale: vector containing the number of scales -- [2] / [2,3] -- if
%     vector size is bigger thant one, data will be repeated to fit the
%     same size
%     *.coarsegrainFeature: cell array containing the name of the feature
%     on which coarsegrain will be applied. Must fit the processing names!
%     [OPTIONAL]
%     -
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.coarsegrainType = 'mean';
% cfg.scale = [15, 30];
% cfg.coarsegrainFeature = {'speed'};
% tsvFile = emcCoarsegrain(tsvFile,cfg);
%
% comments
% the output is written in a field named [featureChosen + 'Cooarsegrain']
% 
% see also
% -
%
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% coarsegrainType: str type of computation of the data -- 'mean' / 'sum'.
% The computation is applied after the segmentation.
errorIfNotField(cfg,'coarsegrainType');
% scale: vector containing the number of scales -- [2] / [2,3] 
errorIfNotField(cfg,'scale');
% coarsegrainFeature: cell array containing the name of the feature on
% which coarsegrain will be applied. Must fit the processing names!
errorIfNotField(cfg,'coarsegrainFeature');
if ~iscell(cfg.coarsegrainFeature)
    cfg.coarsegrainFeature = {cfg.coarsegrainFeature};
end
%% COMPUTATION AREA
for featIdx = 1:numel(cfg.coarsegrainFeature)
    feature = cfg.coarsegrainFeature{featIdx};
    % Check
    if ~isfield(tsvFile.processing, feature)
        error(['Feature ', feature, ' is not part of tsvFile.processing']);
    end
    % Check if multiple columns
    if size(tsvFile.processing.(feature),2) > 1
        data = mean(tsvFile.processing.(feature),2);
    else
        data = tsvFile.processing.(feature);
    end
    if numel(cfg.scale) == 1
        scale = cfg.scale;
        dataCoarse = coarsegrain(data,scale,cfg.coarsegrainType);
    else
        dataCoarse = [];
        for i = 1:numel(cfg.scale)
            scale = cfg.scale(i);
            dataCoarseTmp = coarsegrain(data,scale,cfg.coarsegrainType);
            dataCoarseTmp = repmat(dataCoarseTmp,scale,1); % Repeat the values
            dataCoarseTmp = reshape(dataCoarseTmp,1,numel(dataCoarseTmp));
            % Add to the main dataCoarse -- make it fit the dimension
            if ~isempty(dataCoarse)
                minSize = min(size(dataCoarseTmp,2),size(dataCoarse,2));
                dataCoarse = vertcat(dataCoarse(:,1:minSize),dataCoarseTmp(1:minSize));
            else
                dataCoarse = dataCoarseTmp;
            end
        end
    end
    tsvFile.processing.([feature, 'Coarsegrain']) = dataCoarse;
end
