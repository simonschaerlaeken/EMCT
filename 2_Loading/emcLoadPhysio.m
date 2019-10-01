function [ tsvFile ] = emcLoadPhysio(tsvFile, path, cfg )
%EMCLOADPHYSIO Summary of this function goes here
%   Detailed explanation goes here

%% CHECKING AREA
% type - type of physio recorded - "empathica"
errorIfNotField(cfg, 'physioType')

%% COMPUTATION AREA
if strcmp(cfg.physioType, 'empathica'),
    tsvFile = emcLoadEmpathica( tsvFile, path );
elseif strcmp(cfg.physioType, 'emg'),
    tsvFile = emcLoadEMG( tsvFile, path );
end

end

