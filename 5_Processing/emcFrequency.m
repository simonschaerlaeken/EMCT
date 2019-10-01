function [ tsvFile ] = emcFrequency( tsvFile, cfg )
% Computes a frequency analysis on the movement
% 
% syntax
% tsvFile = emcFrequency(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.deriv: int order of the derivative to be calculated - 1 (speed), 2
%     (acceleration), 3 (jerk)
%     *.euclidianFlag: boolean to determine if the script uses euclidian
%     norm after derivative
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.euclidianFlag = true;
% cfg.display = true
% tsvFile = emcFrequency(tsvFile, cfg)
% 
% comments
% feature is saved in tsvFile.processing.frequency
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end
% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end


%% COMPUTATION AREA
disp('[PROCESSING] Frequency Analysis')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
% Norm
tsvFileNorm = mcnorm(tsvFileMarker);
matrixNorm = tsvFileNorm.data;

Fs = tsvFile.freq;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(matrixNorm,1);             % Length of signal
t = (0:L-1)*T;        % Time vector

n = 2^nextpow2(L);
frequency = fft(matrixNorm,n,1);
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% P1 based on P2 and the even-valued signal length L.
P2 = abs(frequency/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
f = Fs*(0:(L/2))/L;

% Display
if cfg.display
    figure('Name','Frequency Norm','NumberTitle','off')
    imagesc(P1)
end

tsvFile.processing.frequency = f;

end

