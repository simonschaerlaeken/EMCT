function [ ] = emcSaveFeature(varargin)
% Save feature into csv files format (one feature per file) &
% a summary file
% 
% syntax
% emcSaveFeature(tsvFile1[, tsvFile2,...], cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     * cfg.outputDir: str path to output dir
%     * if saveSummaryFeatureFlag = True
%           %     *.summaryType: cell array containing the different type of operation
%                 to be executed while saving the summary file - Features can be either
%                 AVG (averaged), MEDIAN (median), STD (standard deviation), TRIMMEAN
%                 (trim mean), MAX (max value), SUM (summed)
%     [OPTIONAL]
%     *.saveSingleFeatureFlag: boolean flag defining is single feature csv
%     files should be saved
%     *.saveSummaryFeatureFlag: boolean flag defining is summary csv
%     file should be saved
%     *.filename = string containing the prefix name to be used for the
%     summary csv
%     *.saveFeatureList = cell array of all the features to be saved (default: all available)
%     *.saveTrimmeanParam: int value of the cutoff for trimmean (see trimmean)
%     *.summaryPlane: cell array defining the different plane to compute
%     the summary on. Creates different summary for each plane. IN
%     DEVELOPMENT
%     
% output
% -
% 
% examples
% cfg.saveSingleFeatureFlag = true;
% cfg.saveSummaryFeatureFlag = true;
% cfg.saveFeatureList = {'speed', 'convexhull'};
% cfg.summaryType = {'AVG', 'TRIMMEAN'};
% cfg.saveTrimmeanParam = 40;
% emcSaveFeature(tsvFile1,tsvFile2, cfg);
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
if nargin<2
    error('ERROR: No file or no cfg as input')
end
% cfg in last argument
cfg = varargin{end};
% saveSingleFeatureFlag makes the function save all features individually
if ~isfield(cfg, 'saveSingleFeatureFlag')
    cfg.saveSingleFeatureFlag = true;
end
% saveSummaryFeatureFlag makes the function save a summary of all features
if ~isfield(cfg, 'saveSummaryFeatureFlag')
    cfg.saveSingleFeatureFlag = false;
end
% filename is the name of the output for summary files
if ~isfield(cfg, 'filename')
    cfg.filename = 'Summary';
end
% saveFeatureList contains a list of all the features to extract and save
% saveSummaryFeatureFlag makes the function save a summary of all features
undefinedSaveFeatureList = false;
if ~isfield(cfg, 'saveFeatureList')
    undefinedSaveFeatureList = true;
end
% summaryType - list of all the type of concatenation you want to have -
% AVG, STD, MEDIAN,..
if ~isfield(cfg, 'summaryType')
    disp('[WARNING] No summary type defined. Default: AVG.')
    cfg.summaryType = {'AVG'};
end
% summaryType - list of all the type of concatenation you want to have -
% AVG, STD, MEDIAN,..
if ~isfield(cfg, 'summaryPlane')
    disp('[WARNING] No summary plane defined. Default: {}.')
    cfg.summaryPlane = {};
end

% summaryType - list of all the type of concatenation you want to have -
% AVG, STD, MEDIAN,..
if ~isfield(cfg, 'concatenateAll')
    disp('[WARNING] concatenateAll set to true. Default: {}.')
    cfg.concatenateAll = true;
end

% saveTrimmeanParam - value for trimmean
if ~isfield(cfg, 'saveTrimmeanParam')
    disp('[WARNING] No saveTrimmeanParam defined. Default: 50.')
    cfg.saveTrimmeanParam = 20;
end

% Check for each arg
allFilename = [];
for fileIdx = 1:nargin-1
    % filename - for title
    if ~isfield(varargin{fileIdx}, 'info')||~isfield(varargin{fileIdx}.info, 'filename')
        filename = strsplit(varargin{fileIdx}.filename, '\');
        filename = filename{end};
        varargin{fileIdx}.info.filename = filename;
    end
    allFilename = [allFilename, ' - ',varargin{fileIdx}.info.filename];
    % Check if "processing" is field
    if ~isfield(varargin{fileIdx}, 'processing')
        error([filename, ' does not have PROCESSING as a field in its structure'])
    else
        % If features are not defined, define them from the first tsvFile
        if undefinedSaveFeatureList
            cfg.saveFeatureList = fieldnames(varargin{fileIdx}.processing);
            undefinedSaveFeatureList = false;
        % Check all the features for all the files
        else
            for featureIdx = 1:numel(cfg.saveFeatureList)
                if ~isfield(varargin{fileIdx}.processing, cfg.saveFeatureList{featureIdx})
                    error([varargin{fileIdx}.info.filename, ' does not have ', cfg.saveFeatureList{featureIdx}, ' as a computed feature.'])
                end
            end
        end
    end
end

%% COMPUTATION AREA
disp(['[SAVE] Saving features for ', allFilename])

% Print one file per feature
if cfg.saveSingleFeatureFlag
    % For every TSV file
    for fileIdx = 1:nargin-1
        tsvFile = varargin{fileIdx};
        for featIdx = 1:length(cfg.saveFeatureList)
            disp(['Writing ', cfg.saveFeatureList{featIdx}])
            % Define Name for file
            featureName = cfg.saveFeatureList{featIdx};
            fileName = [cfg.outputDir, filesep, tsvFile.info.filename '_' featureName '.csv'];
            % Assign
            feature = tsvFile.processing.(featureName);
            % Write Table
            csvwrite(fileName,feature);
        end
    end
end


% Print summary
% Saving decisions are used to defined which feature containing
% multiple columns are going to be contracted into one column (averaged over markers, dimensions...)
if cfg.saveSummaryFeatureFlag
    featureLine = {};
    if cfg.concatenateAll
        savingDecisionsContracted = ones(1, length(cfg.saveFeatureList));
        featureLine = cfg.saveFeatureList;
    else
        savingDecisionsContracted = zeros(1, length(cfg.saveFeatureList));
        for featIdx = 1:numel(cfg.saveFeatureList)
            featureName = cfg.saveFeatureList{featIdx};
            % Check if multiple columns
            if size(varargin{1}.processing.(featureName),2) > 1 % More than 1 column
    %             resp = questdlg({['Should the feature "' cfg.saveFeatureList{featIdx} '" be contracted?']},...
    %                         'Feature Contraction','NO','YES','YES');
                resp = 'YES';
                % If YES, then just written in one column
                if strcmp(resp, 'YES')
                    savingDecisionsContracted(featIdx) = 1;
                    featureLine{end+1} = featureName;
                % If NO, need multiple columns
                else
                    for i = 1:size(eval(['tsvFile.processing.',featureName]),2) %% TODO -- FIX THIS WHEN MULTIPLE
                        featureLine{end+1} = featureName;
                    end
                end
            else % If Feature contains only 1 columns
                featureLine{end+1} = featureName;
            end
        end
    end
    
    % Go through all the summary types
    for summaryTypeIdx = 1:length(cfg.summaryType)
        % Define if we use or not the decomposition by plane
        if ~isempty(cfg.summaryPlane)
            lengthPlane = length(cfg.summaryPlane);
            planeFlag = true;
        else
            lengthPlane = 1;
            planeFlag = false;
        end
        % For each plane (- if not define, just 1 plane which covers every
        % marker)
        for planeIdx = 1:lengthPlane
            % Add suffix or not
            if planeFlag
                planeName = ['-' cfg.summaryPlaneName{planeIdx}];
                planeIdx = findIndexList(varargin{1}.markerName, cfg.summaryPlane{planeIdx}); % a group of the markers
            else
                planeName = '';
                planeIdx = 1:length(varargin{1}.markerName); % all of them
            end
            % Define summary type
            summaryType = cfg.summaryType{summaryTypeIdx};
            disp(['Writing Summary ', summaryType, planeName])
            % Define Name for file
            fileName = [cfg.outputDir, filesep, cfg.filename, '_', summaryType, planeName, '.csv'];
            % Open File
            fid = fopen(fileName , 'w') ;
            if fid == -1
                warning(['File ', fileName, ' might be already opened'])
            end
            % Write the first line
            fprintf(fid, 'Filename,');
            % Write feature line
            for featureNameIdx = 1:length(featureLine)
                fprintf(fid,'%s,',featureLine{featureNameIdx});
            end
            
            % Last columns to be printed are linked to classification
            classificationNameList = fieldnames(varargin{1}.info.classification);
            for classificationIdx = 1:numel(classificationNameList)
                fprintf(fid, '%s,', classificationNameList{classificationIdx});
            end
            fprintf(fid, '\n');
            
            % Write summary feature
            for fileIdx = 1:nargin-1
                tsvFile = varargin{fileIdx};
                fprintf(fid, '%s,', tsvFile.info.filename);
                for featIdx = 1:numel(cfg.saveFeatureList)
                    featureName = cfg.saveFeatureList{featIdx};
                    if length(tsvFile.processing.(featureName)) == 1 % Integer or float but no vector
                        value = tsvFile.processing.(featureName);
                    else
                        % Define the marker to be selected
                        if size(tsvFile.processing.(featureName),2) == 1 % Only one column, can't use the plane
                            markerSelect = 1;
                        else
                            markerSelect = planeIdx;
                            if numel(planeIdx) ~= size(tsvFile.processing.(featureName),2)%% WATCH OUT FOR NOT EQUAL NUM MARKER FOR SEGMENT
                                markerSelect = 1:size(tsvFile.processing.(featureName),2);
                            end
                        end
                        if strcmp(featureName, 'cumuldist') % Special case, takes the last value
                            value = tsvFile.processing.(featureName)(end,:); % TODO - PlaneIdx is not used anymore
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                % According to the mean of contraction
                                if strcmp(summaryType, 'AVG')
                                    value = nanmean(value);
                                elseif strcmp(summaryType, 'MEDIAN')
                                    value = nanmedian(value);
                                elseif strcmp(summaryType, 'TRIMMEAN')
                                    value = trimmean(value, cfg.saveTrimmeanParam);
                                elseif strcmp(summaryType, 'STD')
                                    value = nanstd(value);
                                elseif strcmp(summaryType, 'MAX')
                                    value = max(value);
                                elseif strcmp(summaryType, 'SUM')
                                    value = nansum(value);
                                elseif strcmp(summaryType, 'MODE')
                                    value = roundsd(value,3);
                                    value = mode(value);
                                end
                            end
                        elseif strcmp(summaryType, 'AVG')
                            value = nanmean(tsvFile.processing.(featureName)(:,markerSelect),1);
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = nanmean(value,2);
                            end
                        elseif strcmp(summaryType, 'MEDIAN')
                            value = nanmedian(tsvFile.processing.(featureName)(:,markerSelect),1);
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = nanmedian(value,2);
                            end
                        elseif strcmp(summaryType, 'TRIMMEAN')
                            value = trimmean(tsvFile.processing.(featureName)(:,markerSelect), cfg.saveTrimmeanParam, 1);
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = trimmean(value, cfg.saveTrimmeanParam);
                            end
                        elseif strcmp(summaryType, 'STD')
                            value = nanstd(tsvFile.processing.(featureName)(:,markerSelect));
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = nanmean(value,2);
                            end
                        elseif strcmp(summaryType, 'MAX')
                            value = max(tsvFile.processing.(featureName)(:,markerSelect));
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = max(value,[],2);
                            end
                        elseif strcmp(summaryType, 'SUM')
                            value = nansum(tsvFile.processing.(featureName)(:,markerSelect),1);
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = nansum(value,2);
                            end
                        elseif strcmp(summaryType, 'MODE')
                            value = roundsd(tsvFile.processing.(featureName)(:,markerSelect),3);
                            value = mode(value);
                            if savingDecisionsContracted(featIdx) == 1 % Contracted
                                value = nanmedian(value,2);
                            end
                        end
                    end
                    fprintf(fid, '%d,', value);
                end
                % Find the numerical value for the class referenced
                for classificationIdx = 1:numel(classificationNameList)
                    fprintf(fid, '%s,', tsvFile.info.classification.(classificationNameList{classificationIdx}));
                end
                % The last feature in the list must include the carriage \n
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
    end
end

end

