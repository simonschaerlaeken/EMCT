function [ space, baricenterPoint ] = hmCreate(baricenterPoint, radius, feature)
%HM_CREATE Summary of this function goes here
%   Detailed explanation goes here
% Create empty space from the coordinate of baricenter (Max and min of each axis)
maxX = max(baricenterPoint(:,1));
minX = min(baricenterPoint(:,1));
maxY = max(baricenterPoint(:,2));
minY = min(baricenterPoint(:,2));
width = maxX-minX;
height = maxY-minY;
space = zeros(width,height); % +1 because on Matlabl, matrix starts with index 1
% Shift baricenter so that the minimum is at 1
baricenterPoint(:,1) = baricenterPoint(:,1) - minX;
baricenterPoint(:,2) = baricenterPoint(:,2) - minY;
% fill in space
for t = 1:size(baricenterPoint,1)
    for rangeX = baricenterPoint(t,1)-radius:baricenterPoint(t,1)+radius
        if rangeX > 0 && rangeX <= width
            for rangeY = baricenterPoint(t,2)-radius:baricenterPoint(t,2)+radius
                if rangeY > 0 && rangeY <= height
                    normVect = norm([baricenterPoint(t,1)-rangeX,baricenterPoint(t,2)-rangeY]);
                    if normVect < radius
                        if ~isempty(feature)
                            coef = mean(feature(t,:));
                        else
                            coef = 1;
                        end
                        space(rangeX,rangeY) = space(rangeX,rangeY) + 1*coef/(normVect+1); % modulated by feature
                    end
                end
            end
        end
    end
end

end

