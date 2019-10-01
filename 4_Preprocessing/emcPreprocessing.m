function [ tsvFile ] = emcPreprocessing(tsvFile, cfg)
% Preprocessing of a tsvFile. Includes methods such as mc2frontal, mcrotate, mccenter, mccenter_marker, mcnorm, and mcsmoothen
% from the Motion Capture Toolbox, Copyright 2008, University of Jyvaskyla, Finland 
% 
% syntax
% tsvFile = emcPreprocessing(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY] 
%     *.preprocessingType: cell array with the type(s) of preprocessing 
%       method to be applied
%     if cfg.preprocessingType: 'rotation frontal'
%         *.markerFrontalList: list with the markers index defining the frontal plane (see mc2frontal)
%     if cfg.preprocessingType: 'rotation'
%         *.rotationAngleValue: value of the rotation (see mcrotate)
%         *.rotationAngleAxis: list defining the axis (see mcrotate)
%     if cfg.preprocessingType: 'center marker'
%         *.markerToCenter: string marker name of the marker on which the rest of the data is centered
%     if cfg.preprocessingType: 'smoothen'
%         *.smoothValue: int smooth value (see mcsmoothen)
%     [OPTIONAL] 
%     -
% 
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.preprocessingType = {'rotation frontal', 'rotation', 'center', 'center marker', 'normalize', 'smoothen'};
% cfg.markerFrontalList = [1 2];
% cfg.rotationAngleValue = 30;
% cfg.rotationAngleAxis = [1 0 0];
% cfg.markerToCenter = 'head';
% cfg.smoothValue = 20;
% tsvFile = emcPreprocessing(tsvFile, cfg);
% 
% comments
% Stores the absolute position in tsvFile.info.absolutePosition
% 
% see also
% mc2frontal
% mcrotate
% mccenter
% mccenter_marker
% mcnorm
% mcsmoothen
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
% fillGapFlag - If TRUE, will fill the gap in the motion capture file by
% extrapoling values before and after the gap
if ~isfield(cfg,'preprocessingType')||isempty(cfg.preprocessingType)
    error('ERROR: cfg.preprocessingType not set up. Please set up cfg.preprocessing.')
end
% classificationFlag - if TRUE, classification will be made 
if ~iscell(cfg.preprocessingType)&&isstring(cfg.preprocessingType)
    disp('INFO: cfg.preprocessing is advised to be in a cell shape. Formating done!')
    cfg.preprocessingType = {cfg.preprocessingType};
end

%% COMPUTATION AREA
% Preprocess
for i = 1:length(cfg.preprocessingType)
    % For every preprocessing step
    preprocessing = cfg.preprocessingType{i};
    
    % rotation frontal - rotate the motion capture data along an axis drawn
    % with two markers defined in cfg.markerFrontalList
    if strcmp('rotation frontal', preprocessing)
        % CHECKING AREA
        errorIfNotField(cfg,'markerFrontalList')
        % COMPUTATION AREA
        disp('[PREPROCESSING] rotation frontal')
        % Frontal plane: marker 1 & 2
        markerFrontalIdx = findIndexList(tsvFile.markerName,cfg.markerFrontalList);
        % Rotate on frontal plane
        tsvFile = mc2frontal(tsvFile, markerFrontalIdx(1), markerFrontalIdx(2));
        
    % rotation - rotate the motion capture data in a desired angle defined
    % in cfg.rotationAngleValue
    elseif strcmp('rotation', preprocessing)
        % CHECKING AREA
        % Value
        errorIfNotField(cfg,'rotationAngleValue')
        if ~isnumeric(cfg.rotationAngleValue)
            error('ERROR: The angle Axis for the rotation must be numeric')
        end
        % Axis
        errorIfNotField(cfg,'rotationAngleAxis')
        if ~isvector(cfg.rotationAngleAxis)
            error('ERROR: The angle Axis for the rotation must be in a vector format')
        end
        % COMPUTATION AREA
        disp('[PREPROCESSING] rotation')
        for angle_idx = 1:numel(cfg.rotationAngleValue)
            tsvFile = mcrotate(tsvFile, cfg.rotationAngleValue(angle_idx), cfg.rotationAngleAxis(angle_idx,:));
        end
        
    % center - center the data around [0,0,0]    
    elseif strcmp('center', preprocessing)
        disp('[PREPROCESSING] center')
        tsvFile = mccenter(tsvFile);
        
    % center marker - center the data around a specific marker and place that marker at [0,0,0]      
    elseif strcmp('center marker', preprocessing)
        % CHECKING AREA
        % Value
        errorIfNotField(cfg,'markerToCenter')
        % COMPUTATION AREA
        disp('[PREPROCESSING] rotation frontal')
        markerIdx = find(ismember(tsvFile.markerName,cfg.markerToCenter));
        tsvFile = mccenter_marker(tsvFile, markerIdx);
    
    % normalize - normalize the motion capture data
    elseif strcmp('normalize', preprocessing)
        disp('[PREPROCESSING] normalize')
        tsvFile = mcnorm(tsvFile);
        
    % smoothen - smoothen the motion capture data
    elseif strcmp('smoothen', preprocessing)
        % CHECKING AREA
        % Value
        errorIfNotField(cfg,'smoothValue')
        % COMPUTATION AREA
        disp('[PREPROCESSING] smoothen')
        tsvFile = mcsmoothen(tsvFile,cfg.smoothValue);
        

%     elseif strcmp('scrambled', preprocessing) %% NEEEED REWORK!!
%         % Change names of markers
%         tsvFile = scrambled_struct(tsvFile,cfg_main.prefix_body_markers{body_idx});
%         % Rotate
%         tsvFile = mcrotate(tsvFile, 180, [1 0 0]);
%         % ELSE
    else
        disp(['No preprocessing step called ' preprocessing ' exits. Please enter valid preprocessing steps'])
    end
end


end

