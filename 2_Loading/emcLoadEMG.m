function [ tsvFile ] = emcLoadEMG(tsvFile, cfg )
%EMCLOADEMG Function to load the data from EMG files into the TSV;
%Create a field physioData and physioMarkerName

%% CHECKING AREA
% path - refers to the path (directory of file) 
errorIfNotField(cfg,'path')

% Location check
if ~exist(cfg.path)
    error('File does not exist')
elseif isdir(cfg.path) % Directory
    error('emcLoadEMG does not support folder yet')
end

% pairElectrode - pair of electorde to create differential between the two 
errorIfNotField(cfg,'pairElectrode')
if size(cfg.pairElectrode,2) ~= 2, error('Matrix size does not fit pairs'); end
% pairElectrodeName - pair of electorde to create differential between the two 
errorIfNotField(cfg,'pairElectrodeName')
if numel(cfg.pairElectrodeName) ~= size(cfg.pairElectrode,1), error('Names does not have same size as pairs'); end

% epochFlag - says if epoch or not
if ~isfield(cfg,'epochFlag')
    cfg.epochFlag = false;
end

%% COMPUTATION AREA
% Load
emgData = pop_biosig(cfg.path, reshape(cfg.pairElectrode,[1, numel(cfg.pairElectrode)]));

emgData = eeg_checkset( emgData );

% Create electrode differential
emgDataDiff = [];
for i = 1:size(cfg.pairElectrode,1)
    diff = emgData.data(cfg.pairElectrode(i,2),:) - emgData.data(cfg.pairElectrode(i,1),:);
    emgDataDiff = vertcat(emgDataDiff, diff);
end
emgData.data = emgDataDiff;

%Filter data
emgData = pop_eegfiltnew(emgData, 20, 400, 1352, 0, [], 0); % systematiser calcul coefficients filtres
emgData = eeg_checkset( emgData );
% Absolute value
emgData.data = abs(emgData.data);
% Keep only envelope
emgData = pop_eegfiltnew(emgData, [], 40, 676, 0, [], 0);


% save matrix of triggers (from original data EEG)
trigList=[emgData.urevent(:,:).type];
trigLatenc=[EEG.urevent(:,:).latency];
trigger=[trigList; trigLatenc];
trigger=trig';

% Save it
tsvFile.physioData = emgData.data';
tsvFile.physioMarkerName = cfg.pairElectrodeName;
tsvFile.info.physioFreq = 2024;
tsvFile.info.trigger = trigger;

% Epoch
if cfg.epochFlag
    tsvFile.physioData = emcEpoch(tsvFile, cfg);
end
end
