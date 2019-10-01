function [ indexList ] = findIndexListStartsWith(listOriginal, prefix)
%FIND_INDEX_LIST Find the index of the element in a list
% of strings starting with a certain prefix
% [ index_list ] = find_index_list(list_original, list_target)
indexList = [];
for indexListOriginal = 1:numel(listOriginal)
    % Find is string start with prefix
    if startsWith(listOriginal{indexListOriginal},prefix)
        indexList(end+1) = indexListOriginal;
    end
end
if isempty(indexList)
    error(['Prefix ' prefix ' is not part of the marker list. Check for spelling mistakes.']);
end
end
