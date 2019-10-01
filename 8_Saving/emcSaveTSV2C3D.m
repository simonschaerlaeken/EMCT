function emcSaveTSV2C3D( tsvFile, filename )
% Transforms a TSV file into a C3D file ready for export (and export it)
% 
% syntax
% emcSaveTSV2C3D( tsvFile, filename );
% 
% input parameters
% tsvFile: MoCap data structure
% filename: str file name
%     
% output
% 
% 
% examples
% emcSaveTSV2C3D( tsvFile, 'c3dFile.c3d' );
% 
% comments
% -
% 
% see also
% writeC3D.m
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%%
% Change fill the NaN:
if sum(sum(isnan(tsvFile.data)))>0
    tsvFile = mcfillgaps(tsvFile);
end

% Init
c3dFile = struct;

% Markers
marke_data_tsv = tsvFile.data;
Markers = zeros(size(marke_data_tsv,1),size(marke_data_tsv,2)/3,3);
for i = 1:size(marke_data_tsv,2)
    m = ceil(i/3);
    n = mod(i-1,3)+1;
    col = marke_data_tsv(:,i);
    Markers(:,m,n) = col;
end
c3dFile.Markers = Markers;

% VideoFrameRate
VideoFrameRate = tsvFile.freq;
c3dFile.VideoFrameRate = VideoFrameRate;

% Event
if isfield(tsvFile.other,'event')
    Event = tsvFile.other.event;
else
    Event = [];
end
c3dFile.Event = Event;

% ResidualError
if isfield(tsvFile.other,'residualerror')
    ResidualError = repmat(tsvFile.other.residualerror,1,tsvFile.nMarkers);
else
    ResidualError = zeros(size(marke_data_tsv,1),tsvFile.nMarkers);
end
c3dFile.ResidualError = ResidualError;

% AnalogSignals
if isfield(tsvFile,'analogdata')
    AnalogSignals = tsvFile.analogdata;
else
    AnalogSignals = [];
end
c3dFile.AnalogSignals = AnalogSignals;

% AnalogFrameRate
if isempty(AnalogSignals)
    AnalogFrameRate = 0;
elseif isfield(tsvFile.other,'parametergroup')
    AnalogFrameRate = tsvFile.other.parametergroup.ANALOG.RATE.data;
else
    AnalogFrameRate = VideoFrameRate;
end
c3dFile.AnalogFrameRate = AnalogFrameRate;

% CameraInfo
if isfield(tsvFile.other,'camerainfo')
    if isempty(tsvFile.other.camerainfo)
        CameraInfo = zeros(size(ResidualError));
    else
        CameraInfo = tsvFile.other.camerainfo;
    end
else
    CameraInfo = zeros(size(ResidualError));
end
c3dFile.CameraInfo = CameraInfo;

% ParameterGroup
ParameterGroup = struct;
to_delete = {'id','islock','description'};
% Group
if isfield(tsvFile.other,'parametergroup')
    field_name_group = fieldnames(tsvFile.other.parametergroup);
    for group_idx = 1:numel(field_name_group)
        ParameterGroup(group_idx).name = {field_name_group{group_idx}};
        ParameterGroup(group_idx).description = {tsvFile.other.parametergroup.(field_name_group{group_idx}).description};
        % Parameter
        ParameterGroup(group_idx).Parameter = struct;
        field_name_param = fieldnames(tsvFile.other.parametergroup.(field_name_group{group_idx}));
        % Delete useless field names - id, is lock, description
        for to_delete_idx = 1:numel(to_delete)
            idx = find(ismember(field_name_param,to_delete{to_delete_idx}));
            field_name_param(idx) = [];
        end
        % For each param
        for param_idx = 1:numel(field_name_param)
            % name
            ParameterGroup(group_idx).Parameter(param_idx).name = {field_name_param{param_idx}};
            % datatype
            ParameterGroup(group_idx).Parameter(param_idx).datatype = tsvFile.other.parametergroup.(field_name_group{group_idx}).(field_name_param{param_idx}).datatype;
            % datas
            data = tsvFile.other.parametergroup.(field_name_group{group_idx}).(field_name_param{param_idx}).data;
            if ischar(data)
                data = chartable2cellarray(data);
            end
            ParameterGroup(group_idx).Parameter(param_idx).data = data;
            % description
            ParameterGroup(group_idx).Parameter(param_idx).description = {tsvFile.other.parametergroup.(field_name_group{group_idx}).(field_name_param{param_idx}).description};
            % dim
            ParameterGroup(group_idx).Parameter(param_idx).dim = tsvFile.other.parametergroup.(field_name_group{group_idx}).(field_name_param{param_idx}).DIMsize';
        end
    end
else
    ParameterGroup = create_ParameterGroup(tsvFile);
end
c3dFile.ParameterGroup = ParameterGroup;

% -----------------------------------------
% % Export
byteswritten = writeC3D(Markers,VideoFrameRate, ...
				 AnalogSignals, AnalogFrameRate, ...
				 Event,ParameterGroup,CameraInfo, ...
				 ResidualError, filename);


function cellarray = chartable2cellarray(chartable)
cellarray = cell(1,size(chartable,2));
for cellarray_idx = 1:size(chartable,2)
    cell_content = '';
    for j = 1:size(chartable,1)
        char_tmp = chartable(j,cellarray_idx);
        if ~strcmp(char_tmp,' ')
            cell_content = strcat(cell_content, char_tmp);
        end
        cellarray{cellarray_idx} = cell_content;
    end
end

function ParameterGroup = create_ParameterGroup(tsv_file)
ParameterGroup = struct;
% -------------------------- POINT ---------------------------------------
ParameterGroup(1).name = {'POINT'};
ParameterGroup(1).description = {'3-D point parameters'};
ParameterGroup(1).Parameter = struct;
% Parameter in POINT
% - USED
ParameterGroup(1).Parameter(1).name = {'USED'};
ParameterGroup(1).Parameter(1).datatype = 2;
ParameterGroup(1).Parameter(1).data = tsv_file.nMarkers;
ParameterGroup(1).Parameter(1).description = {'Number of trajectorys'};
ParameterGroup(1).Parameter(1).dim = [];
% - SCALE
ParameterGroup(1).Parameter(2).name = {'SCALE'};
ParameterGroup(1).Parameter(2).datatype = 4;
ParameterGroup(1).Parameter(2).data = -0.0510;
ParameterGroup(1).Parameter(2).description = {'Scaling factor'};
ParameterGroup(1).Parameter(2).dim = [];
% - RATE
ParameterGroup(1).Parameter(3).name = {'RATE'};
ParameterGroup(1).Parameter(3).datatype = 4;
ParameterGroup(1).Parameter(3).data = tsv_file.freq;
ParameterGroup(1).Parameter(3).description = {'3D data frame rate'};
ParameterGroup(1).Parameter(3).dim = [];
% - DATA_START
ParameterGroup(1).Parameter(4).name = {'DATA_START'};
ParameterGroup(1).Parameter(4).datatype = 2;
ParameterGroup(1).Parameter(4).data = 8;
ParameterGroup(1).Parameter(4).description = {'Pointer to first block 3D/analog data'};
ParameterGroup(1).Parameter(4).dim = [];
% - FRAMES
ParameterGroup(1).Parameter(5).name = {'FRAMES'};
ParameterGroup(1).Parameter(5).datatype = 2;
data = size(tsv_file.data,1);
if data > 65535, % Limit because of the datatype (see Long_Frame for longer)
    data = 65535;
end
ParameterGroup(1).Parameter(5).data = data;
ParameterGroup(1).Parameter(5).description = {'Numbers of frames'};
ParameterGroup(1).Parameter(5).dim = [];
% - LABELS
ParameterGroup(1).Parameter(6).name = {'LABELS'};
ParameterGroup(1).Parameter(6).datatype = -1;
ParameterGroup(1).Parameter(6).data = tsv_file.markerName';
ParameterGroup(1).Parameter(6).description = {'Trajectories labels'};
ParameterGroup(1).Parameter(6).dim = [32,numel(tsv_file.markerName)];
% - DESCRIPTIONS
ParameterGroup(1).Parameter(7).name = {'DESCRIPTIONS'};
ParameterGroup(1).Parameter(7).datatype = -1;
ParameterGroup(1).Parameter(7).data = tsv_file.markerName';
ParameterGroup(1).Parameter(7).description = {'Point Descriptions'};
ParameterGroup(1).Parameter(7).dim = [50,numel(tsv_file.markerName)];
% - UNITS
ParameterGroup(1).Parameter(8).name = {'UNITS'};
ParameterGroup(1).Parameter(8).datatype = -1;
ParameterGroup(1).Parameter(8).data = {'mm'};
ParameterGroup(1).Parameter(8).description = {'Measurment units'};
ParameterGroup(1).Parameter(8).dim = 2;
% - LONG_FRAMES
ParameterGroup(1).Parameter(9).name = {'LONG_FRAMES'};
ParameterGroup(1).Parameter(9).datatype = 4;
ParameterGroup(1).Parameter(9).data = size(tsv_file.data,1);
ParameterGroup(1).Parameter(9).description = {'Number of frames (not limited to 65535)'};
ParameterGroup(1).Parameter(9).dim = [];
% - Y_SCREEN
ParameterGroup(1).Parameter(10).name = {'Y_SCREEN'};
ParameterGroup(1).Parameter(10).datatype = -1;
ParameterGroup(1).Parameter(10).data = {'-Z'};
ParameterGroup(1).Parameter(10).description = {'Axis of screen Y coordinate'};
ParameterGroup(1).Parameter(10).dim = 2;
% - X_SCREEN
ParameterGroup(1).Parameter(11).name = {'X_SCREEN'};
ParameterGroup(1).Parameter(11).datatype = -1;
ParameterGroup(1).Parameter(11).data = {'-X'};
ParameterGroup(1).Parameter(11).description = {'Axis of screen X coordinate'};
ParameterGroup(1).Parameter(11).dim = 2;

% -------------------------- ANALOG ---------------------------------------
ParameterGroup(2).name = {'ANALOG'};
ParameterGroup(2).description = {'Analog data parameters'};
ParameterGroup(2).Parameter = struct;
% Parameter in ANALOG
% - USED
ParameterGroup(2).Parameter(1).name = {'USED'};
ParameterGroup(2).Parameter(1).datatype = 2;
if isfield(tsv_file,'analogdata')
    data = size(tsv_file.analogdata,2);
else
    data = 0;
end
ParameterGroup(2).Parameter(1).data = data;
ParameterGroup(2).Parameter(1).description = {'Number of analog channels'};
ParameterGroup(2).Parameter(1).dim = [];
% - LABELS
ParameterGroup(2).Parameter(2).name = {'LABELS'};
ParameterGroup(2).Parameter(2).datatype = -1;
if isfield(tsv_file,'analogName')
    data = tsv_file.analogName;
else
    data = [];
end
ParameterGroup(2).Parameter(2).data = data;
ParameterGroup(2).Parameter(2).description = {'Analog labels'};
ParameterGroup(2).Parameter(2).dim = [32,numel(data)];
% - DESCRIPTIONS
ParameterGroup(2).Parameter(3).name = {'DESCRIPTIONS'};
ParameterGroup(2).Parameter(3).datatype = -1;
ParameterGroup(2).Parameter(3).data = data; % Same as for LABELS
ParameterGroup(2).Parameter(3).description = {'Analog Descriptions'};
ParameterGroup(2).Parameter(3).dim = [50,numel(data)];
% - GEN_SCALE
ParameterGroup(2).Parameter(4).name = {'GEN_SCALE'};
ParameterGroup(2).Parameter(4).datatype = 4;
ParameterGroup(2).Parameter(4).data = 1;
ParameterGroup(2).Parameter(4).description = {'Universal analog scaling factor'};
ParameterGroup(2).Parameter(4).dim = [];
% - SCALE
ParameterGroup(2).Parameter(5).name = {'SCALE'};
ParameterGroup(2).Parameter(5).datatype = 4;
ParameterGroup(2).Parameter(5).data = [];
ParameterGroup(2).Parameter(5).description = {'Analog scaling parameters'};
ParameterGroup(2).Parameter(5).dim = 0;
% - OFFSET
ParameterGroup(2).Parameter(6).name = {'OFFSET'};
ParameterGroup(2).Parameter(6).datatype = 2;
ParameterGroup(2).Parameter(6).data = [];
ParameterGroup(2).Parameter(6).description = {'Analog offsets'};
ParameterGroup(2).Parameter(6).dim = 0;
% - UNITS
ParameterGroup(2).Parameter(7).name = {'UNITS'};
ParameterGroup(2).Parameter(7).datatype = -1;
ParameterGroup(2).Parameter(7).data = [];
ParameterGroup(2).Parameter(7).description = {'Analog Measurement units'};
ParameterGroup(2).Parameter(7).dim = [4,0];
% - RATE
ParameterGroup(2).Parameter(8).name = {'RATE'};
ParameterGroup(2).Parameter(8).datatype = 4;
if isfield(tsv_file,'anaFreq')
    data = tsv_file.anaFreq;
else
    data = [];
end
ParameterGroup(2).Parameter(8).data = data;
ParameterGroup(2).Parameter(8).description = {'Analog data sample rate'};
ParameterGroup(2).Parameter(8).dim = [];

% -------------------------- SEG ---------------------------------------
ParameterGroup(3).name = {'SEG'};
ParameterGroup(3).description = {'Seg data parameters'};
ParameterGroup(3).Parameter = struct;
% Parameter in POINT
% - PRED_ERROR
ParameterGroup(3).Parameter(1).name = {'PRED_ERROR'};
ParameterGroup(3).Parameter(1).datatype = 4;
ParameterGroup(3).Parameter(1).data = 30;
ParameterGroup(3).Parameter(1).description = {'Maximum predictor error'};
ParameterGroup(3).Parameter(1).dim = [];
% - MAX_RESID
ParameterGroup(3).Parameter(2).name = {'MAX_RESID'};
ParameterGroup(3).Parameter(2).datatype = 4;
if isfield(tsv_file.other,'residualerror')
    data = max(tsv_file.other.residualerror);
else
    data = 0;
end
ParameterGroup(3).Parameter(2).data = data;
ParameterGroup(3).Parameter(2).description = {'Max residuals'};
ParameterGroup(3).Parameter(2).dim = [];
% - ACC_FACTOR
ParameterGroup(3).Parameter(3).name = {'ACC_FACTOR'};
ParameterGroup(3).Parameter(3).datatype = 4;
ParameterGroup(3).Parameter(3).data = 50000;
ParameterGroup(3).Parameter(3).description = {'Segment starting acceleration'};
ParameterGroup(3).Parameter(3).dim = [];
% - NOISE_FACTOR
ParameterGroup(3).Parameter(4).name = {'NOISE_FACTOR'};
ParameterGroup(3).Parameter(4).datatype = 4;
ParameterGroup(3).Parameter(4).data = 10;
ParameterGroup(3).Parameter(4).description = {'Segment starting noise factor'};
ParameterGroup(3).Parameter(4).dim = [];
% - DATA_LIMITS
ParameterGroup(3).Parameter(5).name = {'DATA_LIMITS'};
ParameterGroup(3).Parameter(5).datatype = 4;
data = [min(min(tsv_file.data(:,1:3:end))),min(min(tsv_file.data(:,2:3:end))),min(min(tsv_file.data(:,3:3:end)));...
        max(max(tsv_file.data(:,1:3:end))),max(max(tsv_file.data(:,2:3:end))),max(max(tsv_file.data(:,3:3:end)))];
ParameterGroup(3).Parameter(5).data = data;
ParameterGroup(3).Parameter(5).description = {'Measurement volume limits'};
ParameterGroup(3).Parameter(5).dim = [2,3];

% -------------------------- MANUFACTURER ---------------------------------------
ParameterGroup(4).name = {'MANUFACTURER'};
ParameterGroup(4).description = {'Manufacturer information'};
ParameterGroup(4).Parameter = struct;
% Parameter in POINT
% - COMPANY
ParameterGroup(4).Parameter(1).name = {'COMPANY'};
ParameterGroup(4).Parameter(1).datatype = -1;
ParameterGroup(4).Parameter(1).dim = 8;
ParameterGroup(4).Parameter(1).data = {'Qualisys'};
ParameterGroup(4).Parameter(1).description = {'Company name'};

% - SOFTWARE
ParameterGroup(4).Parameter(2).name = {'SOFTWARE'};
ParameterGroup(4).Parameter(2).datatype = -1;
ParameterGroup(4).Parameter(2).dim = 22;
ParameterGroup(4).Parameter(2).data = {'Qualisys Track Manager'};
ParameterGroup(4).Parameter(2).description = {'Software used to create measurement file'};

% - VERSION
ParameterGroup(4).Parameter(3).name = {'VERSION'};
ParameterGroup(4).Parameter(3).datatype = 2;
ParameterGroup(4).Parameter(3).dim = 3;
ParameterGroup(4).Parameter(3).data = [2;6;680];
ParameterGroup(4).Parameter(3).description = {'Software version used to create measurement file'};

