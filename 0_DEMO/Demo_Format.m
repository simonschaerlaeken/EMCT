%% Paths
pathComputer = 'C:\Users\schaerla\Google_Drive\PhD'; %Windows
% pathComputer = '/Users/Simon/Google_Drive/PhD'; % MAc

% Add path for toolbox
addpath(genpath([pathComputer, filesep, '0_PROJETS', filesep, '04_Motion_Capture', filesep, '043_Code_Matlab', filesep, '0434_EMC_Toolbox']))
cfg = []; % Empty cfg

%% DIRECTORY
% cfg.inputDir = [pathComputer, filesep, '7_Travail_Collaboratif', filesep, '72_MasterStudents', filesep, '072A_Cyrille', filesep, 'C3D', filesep, 'AddMarker'];
cfg.inputDir = [pathComputer, filesep, '7_Travail_Collaboratif', filesep, '72_MasterStudents', filesep, '072A_Cyrille', filesep, 'C3D'];
cfg.outputDir = [pathComputer, filesep, '7_Travail_Collaboratif', filesep, '72_MasterStudents', filesep, '072A_Cyrille', filesep, 'TSV_Complete'];
%% FILENAMES
% Code to look up all the files in directory and keep those with specific
% extension
directoryInfos = dir(cfg.inputDir);
listFilesC3d = {directoryInfos.name}';
cfg.listFiles = {};
for i = 1:length(listFilesC3d)
    [~,~,ext] = fileparts(listFilesC3d{i});
    if strcmp(ext, '.c3d') % .tsv .c3d
        cfg.listFiles{end+1} = listFilesC3d{i};
    end
end
cfg.listFiles = cfg.listFiles';
%cfg.listFiles = {'BaExAs (1).c3d'};
%% EPOCH DATA
% epochdata = [858,4514;56,3885;473,4303;408,3927;321,3897;545,4153];
epochdata = [463,4127;409,4087;429,3997;560,4074;406,3952;470,3918;707,3945;447,3915;507,3929;680,3818;253,4012;595,3863;355,4192;370,3947;354,3927;355,3903;555,3995;257,3916];
%% LOAD C3D
duration = zeros(numel(cfg.listFiles),1);
for i = 1:numel(cfg.listFiles)
    filename = cfg.listFiles{i};
    disp(filename)
    filename = cfg.listFiles{i};

    cfg.fillGapFlag = false;
    % cfg.deleteMarker = {'*78','*79'};
    %     cfg.removePrefixFlag = false;
    tsvFile = emcLoadC3D([cfg.inputDir, filesep, filename], cfg);
%     cfg.transformMarker = {'Body:LShoulderBack','Body:RShoulderBack'};
%     cfg.transformName = 'Body:BackTop';
%     cfg.betweenParameter = 50;
%     tsvFile = emcMidpoint(tsvFile, cfg);
    %% EPOCH
    %cfg.epochInOut = [epochdata(i,1), epochdata(i,2)];
    %tsvFile = emcEpoch(tsvFile, cfg);
    %     duration(i) = size(tsvFile.data,1)/tsvFile.freq;
    %     %% PROCESS
    %cfg.fillGapFlag = true;
    %cfg.removePrefixFlag = false;
    

    cfg.orderMarkerFlag = true;
    cfg.orderMarker = {'Body:WaistLFront';'Body:WaistRFront';'Body:WaistLBack';'Body:WaistRBack';'Body:BackTop';'Body:Chest';'Body:BackLeft';'Body:BackRight';'Body:HeadTop';'Body:HeadFront';'Body:HeadSide';'Body:LShoulderBack';'Body:LShoulderTop';'Body:LElbowOut';'Body:LUArmHigh';'Body:LHandOut';'Body:LWristOut';'Body:LWristIn';'Body:RShoulderBack';'Body:RShoulderTop';'Body:RElbowOut';'Body:RUArmHigh';'Body:RHandOut';'Body:RWristOut';'Body:RWristIn';'Body:LKneeOut';'Body:LThigh';'Body:LAnkleOut';'Body:LShin';'Body:LToeOut';'Body:LToeIn';'Body:RKneeOut';'Body:RThigh';'Body:RAnkleOut';'Body:RShin';'Body:RToeOut';'Body:RToeIn';'Djembe:Marker1';'Djembe:Marker2';'Djembe:Marker3';'Djembe:Marker4';'Djembe:Marker5'};
    %tsvFile = emcOrderMarker( tsvFile, cfg );
    emcWriteTsv(tsvFile, cfg);
    %% SAVE C3D
    %emcSaveTSV2C3D(tsvFile, [cfg.outputDir,filesep, filename(1:end-3),'c3d'])
end
