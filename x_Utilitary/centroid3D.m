function [centroid_data] = centroid3D(data)
% Function to find the centroid point for n points in 3D
nb_points = size(data,2)/3;
centroid_data = zeros(size(data,1), 3);
for marker_idx = 1:3:size(data,2),
    centroid_data = centroid_data + data(:,marker_idx:marker_idx+2);
end
centroid_data = centroid_data/nb_points;
end
