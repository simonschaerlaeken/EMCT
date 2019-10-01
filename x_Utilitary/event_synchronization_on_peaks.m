function [out_Q, out_q] = event_synchronization_on_peaks(varargin)

% Usage:
%   [out_Q, out_q] = event_synchronization_on_peaks(in_s1, in_s2, ...
%           in_alg, in_window, in_make_plot);
%
% Description:
%   This function computes event synchronization on the peaks of the two 
%   univariate time series it receives as input. Event synchronization 
%   is computed with the algorithm proposed by Quian Quiroga 
%   (Phys.Rev. E 66 041904).
%
% Input arguments:
%   in_s1:          The first input time series. It should be an 
%                   univariate time series provided as a vector.
%   in_s2:          The second input time series. It should be an 
%                   univariate time series provided as a vector. 
%   in_alg:         A string specifying the kind of output to be 
%                   obtained from the Quian Quiroga's algorithm. 
%                   The following options are supported:
%                       - 'tot': the outputs out_Q and out_q measure 
%                          respectively the synchonization (Q(tau)) and 
%                          the delay behavior (q(tau)) between the whole   
%                          time series;
%                       - 'tsl': the outputs out_Q and out_q are the time 
%                          resolved variants of Q(tau) and q(tau);
%                       - 'ssl': the outputs out_Q and out_q are the  
%                          average of Q(tau) and q(tau) over the last  
%                          'window' time steps.
%   in_window       The size of the window (in number of samples) used for
%                   computing the averaged version of Q(tau) and q(tau)
%                   when the 'ssl' option is selected.
%   in_make_plot    If different from zero, a plot of the input time
%                   series and the detected peaks is generated.
%
% Outputs:
%   out_Q:          The synchonization (Q(tau)) between the peaks of the  
%                   input time series. This is a scalar when 'tot' is  
%                   selected or a time series when 'ssl' is selected.
%   out_q:          The delay behavior (q(tau)) between the peaks of the  
%                   input time series. This is a scalar when 'tot' is  
%                   selected or a time series when 'ssl' is selected.
%
% Remarks:
%   in_s1 and in_s2 are the only mandatory parameters.
%   If not specified, the default value for in_alg is 'tot' and the default
%   value for in_window is 75.
%   This function uses the IPEM Toolbox for extraction of peaks and the Qqm
%   function by Giovanna Varni implementing the Quian Quiroga's algorithm.
%
% Example:
%   [Q, q] = event_synchronization_on_peaks(s1, s2, 'ssl', 100, 1);
%
% Authors:
%   Gualtiero Volpe
%   Created on: 20110610
%   Revised on: 20110611

%--------------------------------------------------------------------------
% Copyright (c) 2011 - InfoMus - DIST - University of Genova
%
% http://www.infomus.org
% mailto:info@infomus.org
%
%--------------------------------------------------------------------------
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
% -------------------------------------------------------------------------

available_algs = {'tot', 'tsl', 'ssl'};

[in_s1, in_s2, in_alg, in_window, in_make_plot] = ...
    IPEMHandleInputArguments(varargin, 3, {[], [], 'tot', 75, 0});
% in_s1        = varargin{1};
% in_s2        = varargin{2};
% in_alg       = varargin{3};
% in_window    = 75;
% in_make_plot = varargin{4};

if (isempty(in_s1) || isempty(in_s2))
    error('The input time series are empty. Please check your data.')
end

[m1, n1] = size(in_s1);
if (~((m1 == 1) || (n1 == 1)) || (m1 == 1 && n1 == 1))
    error('The input time series should be provided as a vector.')
end
if (n1 == 1)
    in_s1 = in_s1';
    n1 = size(in_s1, 2);
end

[m2, n2] = size(in_s2);
if (~((m2 == 1) || (n2 == 1)) || (m2 == 1 && n2 == 1))
    error('The input time series should be provided as a vector.')
end
if (n2 == 1)
    in_s2 = in_s2';
    n2 = size(in_s2, 2);
end

peaks1 = IPEMFindAllPeaks(in_s1, [], 0)';
peaks2 = IPEMFindAllPeaks(in_s2, [], 0)';

% If no peaks:
if isempty(peaks1)
    peaks1 = [1];
end
if isempty(peaks2)
    peaks2 = [1];
end

% [peaks1, valleys1] = siPeakDet(in_s1, 0.4);
% [peaks2, valleys2] = siPeakDet(in_s2, 0.4);

%time1 = [peaks1, ones(size(peaks1, 1), 1)];
%time2 = [peaks2, ones(size(peaks2, 1), 1)];

time1 = [peaks1(:, 1), ones(size(peaks1, 1), 1)];
time2 = [peaks2(:, 1), ones(size(peaks2, 1), 1)];

switch (in_alg)
    case available_algs{1}
        [out_Q, out_q] = Qqm2(time1, time2, size(time1, 1), size(time2, 1), available_algs{1});
    case available_algs{2}
        [out_Q, out_q] = Qqm2(time1, time2, size(time1, 1), size(time2, 1), [], available_algs{2}, max(n1, n2));
    case available_algs{3}
        [out_Q, out_q] = Qqm2(time1, time2, size(time1, 1), size(time2, 1), in_window, available_algs{3}, max(n1, n2));
end

if (in_make_plot ~= 0)
    figure;
    subplot(2, 1, 1)
    plot(1 : n1, in_s1, 'b', peaks1, in_s1(peaks1), 'ro')
    subplot(2, 1, 2)
    plot(1 : n2, in_s2, 'b', peaks2, in_s2(peaks2), 'ro')
end

    
