function saen = sampleEntropy( dim, r, data, tau )
% SAMPEN Sample Entropy
%   calculates the sample entropy of a given time series data

%   SampEn is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Mar 21, 2012
%---------------------------------------------------------------------
%% see Sampen
%https://www.researchgate.net/publication/230640364_On_the_relation_between_correlation_dimension_approximate_and_sample_entropy_parameters_and_a_fast_algorithm_for_their_calculation
if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end

N = length(data);
correl = zeros(1,2);
dataMat = zeros(dim+1,N-dim);
for i = 1:dim+1
    dataMat(i,:) = data(i:N-dim+i-1);
end

for m = dim:dim+1
    count = zeros(1,N-dim);
    tempMat = dataMat(1:m,:);
    
    for i = 1:N-m
        % calculate Chebyshev distance, excluding self-matching case
        dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
        
        % calculate Heaviside function of the distance
        % User can change it to any other function
        % for modified sample entropy (mSampEn) calculation
        D = (dist < r);
        
        count(i) = sum(D)/(N-dim);
    end
    
    correl(m-dim+1) = sum(count)/(N-dim);
end

saen = log(correl(1)/correl(2));
end

% COMMENTS
% The DOI of the paper is 10.1109/MEMB.2009.934629 so you can find it
% easily.
% 
% The first argument, dim, means embedded dimension. That is the length of
% the sequence to be compared, and typical value for that is 2 or 3. The
% second argument, r, is the tolerance for matching. In other words, a pair
% of values whose difference is smaller than r is considered matching.
% Typically 0.2*std is used (std means standard deviation of the whole
% timeseries data), but this paper by Chon suggests that r value should be
% chosen carefully, and one should look for r value that maximizes the
% ApEn.
