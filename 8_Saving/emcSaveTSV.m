function [] = emcSaveTSV( tsvFile, path )
% Save a TSV file structure as a tsv file using mcwritetsv (University of Jyvaskyla, Finland)
% 
% syntax
% c3dFile = emcSaveTSV( tsvFile, filename );
% 
% input parameters
% tsvFile: MoCap data structure
% path: str path folder name
%     
% output
% -
% 
% examples
% emcSaveTSV( tsvFile, 'C:/desktop' );
% 
% comments
% -
% 
% see also
% mcwritetsv.m
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
%% COMPUTATION AREA
mcwritetsv(tsvFile, path)
end

