function [ tsvFileStruct ] = emcLoad(filename, cfg)
%EMCLOAD Function loading a tsv file into matlab.
%   Input:
%       path: data path
%       cfg: options (contains: fillGapFlag)

%% CHECKING AREA
% Location check
if ~exist(filename, 'file')
    error('File does not exist')
end
if isdir(filename) % Directory
    error('ERROR: Folder not yet supported, please import files individually. emcFormat can help you set up your files.')
end
if ~any(regexp(filename,'tsv$')) % if not TSV fole
    error('ERROR: Only TSV files supported. emcFormat can help you set up your files (convertion available from C3D).')
end

%% COMPUTATION
tsvFileStruct = emcLoadTSV(filename, cfg);

end

