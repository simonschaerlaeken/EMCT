function [y] = coarsegrain(X,S,type) 
% e.g., coarsegrain(z,5)

%Input Parameters
%X input signal vector
%S number of scales

%Output Parameters
%y coarsegrained series


if S==1 %exceptional solution where coarse-grained signal is the input signal vector 
 y = X';
 return
end

N = length(X);

%insure that vector is multiple of scales ~= 
if (mod(N,S)~=0)
 X=cat(1,X,zeros(S-mod(N,S),1));% fullfil
end
N = length(X);
if strcmp(type, 'sum')
    y = sum(reshape(X,S,N/S)); % split fullfilled original signal into a new matrix with subvector or length S.then Sum it or average.
elseif strcmp(type, 'mean')
    y = mean(reshape(X,S,N/S));
else
    error('This type is not defined in coarsegrain')
end
end

