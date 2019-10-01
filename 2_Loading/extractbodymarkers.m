function [ bodyMarkersList, bodyMarkersIdxList ] = extractbodymarkers( tsvFile, bodyPrefix )
%EXTRACT_BODY_MARKER Function to extract the list of Markers for one body
% referenced in prefix from a TSV file
bodyMarkersList = {};
bodyMarkersIdxList = [];
i = 1;
allMarkersList = tsvFile.markerName;
for markerIdx = 1:length(allMarkersList),
    % If starts with the prefix, add it
    if strncmpi(allMarkersList{markerIdx}, bodyPrefix, length(bodyPrefix))
        bodyMarkersList(1,i) = allMarkersList(markerIdx);
        bodyMarkersIdxList(end+1) = markerIdx;
        i = i + 1;
    end
end
end

