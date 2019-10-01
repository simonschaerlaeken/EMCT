function [ list, success ] = removePrefixList( list, prefix )
%REMOVEPREFIXLIST Remove a defined prefix in every element of a list (Every
% must contain the prefix for success to be true)
listTmp = list;
success = true; % Show if found
for elementIdx = 1:numel(list)
    if strncmp(list{elementIdx}, prefix, length(prefix))
        listTmp{elementIdx} = list{elementIdx}(length(prefix)+1:end);
    end
end
list = listTmp;
end

