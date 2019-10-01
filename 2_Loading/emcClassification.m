function [ tsvFile ] = emcClassification(tsvFile, cfg)
% Adds tags to the tsv structure for classification. This can be made
% according to different methods: "prefix", uses the filename given and
% extract information at specific location
% 
% syntax
% tsvFile = emcClassification(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.classificationType: define the classification method
%     if cfg.classificationType: 'prefix'
%       *.classInfo: cell array containing the classes name
%       *.classPrefixPosition: cell array containing the indexes of the
%       classification letter in the filename
%       *.classListPrefix: cell array containing the different option
%       available at the positions mentionned in cfg.classPrefixPosition
%       *.classListValue: cell array containing the value of each of the
%       options available in cfg.classListPrefix
%     [OPTIONAL]
%     -
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.classificationType = 'prefix';
% cfg.classInfo = {'className1','className2'};
% cfg.classPrefixPosition = {1,2}; 
% cfg.classListPrefix = {{'1','2','3'},...
%                        {'1','2'}};
% cfg.classListValue = {{'class1_1','class1_2','class1_3'},...
%                       {'class1_2','class2_2'}}; 
% tsvFile = emcClassification(tsvFile, cfg);
% 
% 
% 
% comments
% if store the classification in tsvFile.info.classification
% 
% see also
% emcLoadSingle
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
% classificationType - defines which classification type you want to apply
%       - prefix: looks at the prefix of the filename
errorIfNotField(cfg,'classificationType');
if strcmp(cfg.classificationType, 'prefix') % Check argument for prefix
    errorIfNotField(cfg,'classInfo');
    errorIfNotField(cfg,'classPrefixPosition');
    errorIfNotField(cfg,'classListPrefix');
    errorIfNotField(cfg,'classListValue');
end
%% COMPUTATION AREA
filename = strsplit(tsvFile.filename, filesep);
filename = filename{end};
% Prefix mode
if strcmp(cfg.classificationType, 'prefix')
    for classInfoIdx = 1:length(cfg.classInfo)
        classInfo = cfg.classInfo{classInfoIdx};
        prefix = filename(cfg.classPrefixPosition{classInfoIdx});
        classIdx = findIndexList(cfg.classListPrefix{classInfoIdx},{prefix});
        class = cfg.classListValue{classInfoIdx}{classIdx};
        % Create field in tsvFile
        tsvFile.info.classification.(classInfo) = class;
    end
end

end

