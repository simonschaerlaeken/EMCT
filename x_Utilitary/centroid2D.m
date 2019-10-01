function [centroid_data] = centroid2D(data)
% Function to find the centroid point for n points in 2D
nb_points = size(data,2)/2;
centroid_data = zeros(size(data,1), 2);
for marker_idx = 1:2:size(tsv_file.data,2),
    centroid_data = centroid_data + tsv_file.data(:,marker_idx:marker_idx+2);
end
centroid_data = centroid_data/nb_points;
end

