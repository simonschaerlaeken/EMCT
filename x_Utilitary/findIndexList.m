function [ index_list ] = findIndexList(list_original, list_target)
%FIND_INDEX_LIST Find the index of the element from a target list in a list
% of strings
% [ index_list ] = find_index_list(list_original, list_target)
index_list = zeros(1,length(list_target));
for list_target_idx = 1:length(list_target)
    % Find the index of target in orginal list
    idx = find(ismember(list_original,list_target{list_target_idx}));
%     for i = 1:length(list_original)
%         if strcmp(list_original{i},list_target{list_target_idx})
%             idx = i;
%         end
%     end
    %idx = find(not(cellfun('isempty', strfind(list_original, list_target{list_target_idx}))));
    if isempty(idx)
        error(['String ' list_target{list_target_idx} ' is not part of the cell string. Check for spelling mistakes.']);
    end
    index_list(list_target_idx) = idx;
end
end
