function [ dateTime ] = convertDateTime( date, time )
%CONVERTDATETIME Convert date in format "YYYY-MM-DD" and time in format
%"HH:MM:SS" into a datetime structure
%   Detailed explanation goes here
dateSplit = strsplit(date, '/');
timeSplit = strsplit(time, ':');
dateTime = datetime(str2num(dateSplit{3}),str2num(dateSplit{2}),str2num(dateSplit{1}),str2num(timeSplit{1}),str2num(timeSplit{2}),str2num(timeSplit{3}));
end

