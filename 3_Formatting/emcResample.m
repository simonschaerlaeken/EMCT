function [ tsvFile ] = emcResample( tsvFile, cfg )
% Reorders tsvFile placing markers into a specific order
% 
% syntax
% tsvFile = emcResample( tsvFile, cfg );
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.sampleFreqLCM: integer desired frequency for the resampling -- frequency
% correspond to the least common multiple
%     *.dataType:  cell array data name to resample -- 'mocap', 'physio' (,
% 'analog')
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure (with resampled data)
% 
% examples
% cfg.sampleFreqLCM = 64;
% cfg.dataType = {'mocap','physio'};
% tsvFile = emcResample( tsvFile, cfg );
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
% sampleFreq: integer desired frequency for the resampling -- frequency
% correspond to the least common multiple
errorIfNotField(cfg, 'sampleFreqLCM');

% dataType: cell array data name to resample -- 'mocap', 'physio' (,
% 'analog')
errorIfNotField(cfg, 'dataType');
if ~iscell(cfg.dataType)
    cfg.dataType = {cfg.dataType};
end
%% COMPUTATION AREA
disp(['[RESAMPLE] at ', num2str(cfg.sampleFreqLCM), 'Hz'])
for i = 1:numel(cfg.dataType)
    % MOCAP DATA
    if strcmp(cfg.dataType{i}, 'mocap')
        disp('.... Mocap')
        resampleFactor = round(cfg.sampleFreqLCM/tsvFile.freq);
        data = [];
        % For each column since interp can only handle vectors
        for col = 1:size(tsvFile.data,2)
            data = horzcat(data,interp(tsvFile.data(:,col),resampleFactor));
        end
        tsvFile.data = data;
        tsvFile.freq = tsvFile.freq*resampleFactor;
    % PHYSIO DATA
    elseif strcmp(cfg.dataType{i}, 'physio')
        disp('.... Physio')
        for physioIdx = 1:numel(tsvFile.physioFreq)
            % Multiple columns physio
            if size(tsvFile.physioData{physioIdx},2) > 1
                data = [];
                % For each column since interp can only handle vectors
                for col = 1:size(tsvFile.physioData{physioIdx},2)
                    resampleFactor = round(cfg.sampleFreqLCM/tsvFile.physioFreq{physioIdx}(col));
                    data = horzcat(data,interp(tsvFile.physioData{physioIdx}(:,col),resampleFactor));
                    tsvFile.physioFreq{physioIdx}(col) = tsvFile.physioFreq{physioIdx}(col)*resampleFactor;
                end
                tsvFile.physioData{physioIdx} = data;
            % Single column physio
            else
                data = [];
                resampleFactor = round(cfg.sampleFreqLCM/tsvFile.physioFreq{physioIdx});
                data = interp(tsvFile.physioData{physioIdx},resampleFactor);
                tsvFile.physioData{physioIdx} = data;
                tsvFile.physioFreq{physioIdx} = tsvFile.physioFreq{physioIdx}*resampleFactor;
            end
        end
    end
end
end

