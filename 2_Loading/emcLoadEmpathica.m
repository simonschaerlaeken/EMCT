function [ tsvFile ] = emcLoadEmpathica( tsvFile, path)
% Function to load the data from physio into the motion structure
% 
% syntax
% [tsvFile] = emcLoadEmpathica( tsvFile, path )
% 
% input parameters
% tsvFile: motion capture structure
% path: str containign the folder with all the csv files for each physio
%
% output
% tsvFile: motion capture structure
% 
% examples
% [tsvFile] = emcLoadEmpathica( tsvFile, 'C:/physioDir')
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
% Location check
if ~exist(path)
    error('File/Folder does not exist')
elseif isdir(path) % Directory
    % files list
    files = dir(fullfile(path, '*.csv'));
    % Create filename list out of files
    filenameList = cell(1,numel(files));
    for i = 1:numel(files)
        filenameList{i} = fullfile(path,files(i).name);
    end
else
    % If only a file
    filenameList = {path};
end

%% COMPUTATION AREA
disp('[LOAD] Empathica')
for fileIdx = 1:numel(filenameList)
    % Each filename
    filename = filenameList{fileIdx};
    % Extract name of physio
    physioMarkerName = strsplit(filename, filesep);
    physioMarkerName = physioMarkerName{end};
    physioMarkerName = strsplit(physioMarkerName, '.csv');
    physioMarkerName = physioMarkerName{1};
    % Read data
    data = csvread(filename);
    physioTimeStamp = data(1,:);
    physioFreq = data(2,:);
    physioData = data(3:end,:);
    
    % Assign it in tsv
    tsvFile.analogdata(fileIdx).physioName = physioMarkerName;
    tsvFile.analogdata(fileIdx).data = physioData;
    tsvFile.analogdata(fileIdx).freq = physioFreq;
    tsvFile.analogdata(fileIdx).timeStamp = physioTimeStamp;
end

