function [ listModified ] = strcatRecursive(string, list )
%STRCAT_RECURSIVE Apply strcat recursively on List
listModified = list;
for elIdx = 1:length(list)
    if iscellstr(list(elIdx))
        listModified(elIdx) = strcat(string, list(elIdx));
    else
        listModified{elIdx} = strcatRecursive(string, list{elIdx});
    end
end
end