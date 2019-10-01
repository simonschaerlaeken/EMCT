function d2 = mccenter_marker(d,index_marker)
% Translates motion capture data to place a marker at [0 0 0] across markers and over time.
%
% syntax
% d2 = mccenter(d);
%
% input parameters
% d: MoCap structure or data matrix
%
% output
% d2: MoCap structure or data matrix
%
% comments
% Missing data (NaN's) is ignored when calculating the centroid.
%
% Part of the Motion Capture Toolbox, Copyright 2008, 
% University of Jyvaskyla, Finland

d2=[];

if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    start_idx = index_marker*3-2;
    x = mean(mcmean(d.data(:,start_idx)));
    y = mean(mcmean(d.data(:,start_idx+1)));
    z = mean(mcmean(d.data(:,start_idx+2)));
    d2 = mctranslate(d, [-x -y -z]);
elseif isnumeric(d)
    start_idx = index_marker*3-2;
    x = mean(mcmean(d(:,start_idx)));
    y = mean(mcmean(d(:,start_idx+1)));
    z = mean(mcmean(d(:,start_idx+2)));
    d2 = mctranslate(d, [-x -y -z]);
else disp([10, 'The first input argument has to be a variable with MoCap data structure or a data matrix.', 10]);
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
end



end

