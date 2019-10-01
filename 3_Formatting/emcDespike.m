function [ tsvFile ] = emcDespike( tsvFile )
%EMCDESPIKE Summary of this function goes here
%   Detailed explanation goes here

% Calculate difference of position for each point
data1 = tsvFile.data; 
data1 = data1(1:end-1,:);
data2 = tsvFile.data; 
data2 = data2(2:end,:);
diffMat = data2 - data1;
% Find outliers

for i = 1:size(tsvFile.data,2)
%     figure;
%     plot(tsvFile.data(:,i))
%     hold on
    IDXlist = find(diffMat(:,i)==0);
    % MORE than 50% is constant values
    if numel(IDXlist)/size(tsvFile.data,1)>0.9
        % Smoothing
        tsvFile.data(:,i) = smooth(tsvFile.data(:,i),300);
    else
        % Stravinsky Golay
        tsvFile.data(:,i) = sgolayfilt(tsvFile.data(:,i), 3, 301);
    end
end


