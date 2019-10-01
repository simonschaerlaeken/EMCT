filename = '../DATA/VCool/dataOriginal/01_vio_aud_loi_att_mo1_Mocap.csv';
[tsvFile] = emcLoadCSV( filename, cfg);
tsvFile.data = tsvFile.data*1000; % because otherwise it's tiny if you use MOKKA
cfg.transformMarker = {'Head','Back'};
cfg.transformName = 'Neck'; % Creates the neck (could be more precise, but it's ok)
cfg.betweenParameter = 30;
tsvFile = emcBetweenpoint(tsvFile, cfg);
cfg.preprocessingType = {'rotation'};
cfg.rotationAngleValue = 90; % The data is rotated in a weird way, it should put it straight
cfg.rotationAngleAxis = [1 0 0];
tsvFile = emcPreprocessing(tsvFile, cfg);
cfg.connectionMatrix = [1 7; 7 3; 3 4;7 5; 5 6;7 2];
emcPlotBody3D(tsvFile, cfg);
% emcSaveTSV2C3D(test, 'Yo.c3d') %-- If you wanted to save them
% emcSaveTSV(test, pwd) %-- If you wanted to save them