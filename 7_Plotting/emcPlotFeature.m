function [] = emcPlotFeature( tsvFile, cfg )
% Plots the already computed feature
% 
% syntax
% emcPlotFeature(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.featureList: cell array containing the feature names to be plotted
%     [OPTIONAL]
%     *.featureAverageFlag: boolean to specify if you want the average of
%     the feature to appear on the plot
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% -
% 
% examples
% cfg.featureList = {'speed'};
% cfg.featureAverageFlag = true;
% cfg.displayUnit = 's'
% emcPlotFeature(tsvFile, cfg);
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
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% featureList - cell array containing the feature names to be plotted
errorIfNotField(cfg, 'featureList')
if ~iscell(cfg.featureList)
    cfg.featureList = {cfg.featureList};
end
% Check if feature exists
if ~isfield(tsvFile, 'processing')
    error('Processing must be computed before plotting features')
else
    for featureIdx = 1:numel(cfg.featureList)
        featureName = cfg.featureList{featureIdx};
        if ~isfield(tsvFile.processing, featureName)
            error(['Feature ', featureName, 'does not exist in tsvFile.processing'])
        end
    end
end
% featureAverageFlag: boolean to specify if you want the average of the feature to appear on the plot
if ~isfield(cfg, 'featureAverageFlag')
    cfg.featureAverageFlag = false;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
%% WORKING AREA
for featureIdx = 1:numel(cfg.featureList)
    % Feature Name
    featureName = cfg.featureList{featureIdx};
    % Feature Data
    featureData = tsvFile.processing.(featureName);
    % Plot
    % one value
    if size(featureData,1) == 1
        disp(['Feature ', featureName, ' ='])
        disp(featureData)
    % one column
    elseif size(featureData,2) == 1
        figure('Name',[tsvFile.info.filename, ' - ', featureName],'NumberTitle','off');
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
            plot(time, featureData)
        else
            plot(featureData)
        end
        plot(featureData)
    % more columns
    else
        figure('Name',[tsvFile.info.filename, ' - ', featureName],'NumberTitle','off');
        if cfg.featureAverageFlag
            hold on;
            for i = 1:size(featureData,2)
                if strcmp(cfg.displayUnit, 's')
                    time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                    p1 = plot(time, featureData(:,i));
                else
                    p1 = plot(featureData(:,i));
                end
                p1.Color(4) = 0.35;
            end
            plot(mean(featureData,2),'LineWidth',3, 'Color',[0,0.0,0.6]);
            hold off;
        else
            if strcmp(cfg.displayUnit, 's')
                time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                plot(time, featureData);
            else
                plot(featureData);
            end
        end
        if isfield(cfg, 'featMarker') && numel(cfg.featMarker) == size(featureData,2)
            legend(cfg.featMarker);
        elseif numel(tsvFile.markerName) == size(featureData,2)
            legend(tsvFile.markerName);
        end
        
    end
end

