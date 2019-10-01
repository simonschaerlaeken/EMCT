function [ strVect ] = sortStrVect( strVect )
%SORTSTRVECT sort a vector made out of strings
%   Detailed explanation goes here
strVectNum = zeros(size(strVect));
% Tranform str 2 num
for i = 1:numel(strVect)
    strVectNum(i) = str2num(strVect{i})
end
% Sort
strVectNum = sortrows(strVectNum);
% Tranform num 2 str
for i = 1:numel(strVectNum)
    strVect{i} = num2str(strVectNum(i))
end
