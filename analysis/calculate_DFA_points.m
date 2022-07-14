% [DFAobject,DFA_exp] = nbt_Scaling_DFA(DFAobject, Signal, InfoObject,
% CalcInterval, DFA_Overlap);
% Detrended Fluctuation Analysis - Part of the NBT - toolbox
%% Input parameters
% DFAobject         : A DFAobject
% Signal            : The signal
% InfoObject        : The InfoObject that belongs to your signal.
% CalcInterval(1) 	: minimum time-window size computed (in units of seconds!).
% CalcInterval(2) 	: maximum time-window size computed (in units of seconds!).
% DFA_Overlap		: amount of DFA_Overlap between windows (to obtain higher SNR for the fluctuation estimates).

%% output parameters
% DFAobject     : Return the DFAobject with updated information.
% DFA_exp		: The DFA power-law exponent.
%
% Example:
%
% References:
% % Method references...
%
% Peng et al., Mosaic organization of DNA nucleotides, Phys. rev. E (49), 1685-1688 (1994).
% or for a better description: Peng et al., Quantification of scaling exponents and crossover phenomena
% in nonstationary heartbeat time series, Chaos (5), 82-87 (1995).
%
%
% See also:
%   NBT_DFA,

%------------------------------------------------------------------------------------
% Originally created by Klaus Linkenkaer-Hansen (2001), see NBT website (http://www.nbtwiki.net) for current email address
% Improved code - Simon-Shlomo Poil (2008)
% Imported to NBT format. - Simon-Shlomo Poil (2009)
%------------------------------------------------------------------------------------
% 
% ChangeLog - see version control log at NBT website for details.
%
% Copyright (C) 2001  Klaus Linkenkaer-Hansen  (Neuronal Oscillations and Cognition group, 
% Department of Integrative Neurophysiology, Center for Neurogenomics and Cognitive Research, 
% Neuroscience Campus Amsterdam, VU University Amsterdam)
%
% Part of the Neurophysiological Biomarker Toolbox (NBT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% See Readme.txt for additional copyright information.
% -


function [DFA_x,DFA_y] = nbt_doDFA(Signal, Fs, CalcInterval, DFA_Overlap)

    %% Get or set default parameters...

    res_logbin = 10;

    %******************************************************************************************************************
    %% Begin analysis
    % Find DFA_x

    % Defining window sizes to be log-linearly increasing.
    d1 = floor(log10(CalcInterval(1)*Fs));
    d2 = ceil(log10(CalcInterval(2)*Fs));
    DFA_x_t = round(logspace(d1,d2,(d2-d1)*res_logbin));	% creates vector from 10^d1 to 10^d2 with N log-equidistant points.
    DFA_x = DFA_x_t((CalcInterval(1)*Fs <= DFA_x_t & DFA_x_t <= CalcInterval(2)*Fs));	% Only include log-bins in the time range of interest!
    DFA_x = DFA_x';

    %% The DFA algorithm...
    DFA_y = nan(size(DFA_x,1),1);
    %% Do check of CalcInterval
    %         if (CalcInterval(2) > 0.1*length(Signal(:,GetChannelID))/Fs)
    %             display('The upper limit of CalcInterval is larger than the recommended 10% of the signal lenght')
    %         end

    y = Signal(:,1)./mean(Signal(:,1));
    % First we convert the time series to a series of fluctuations y(i) around the mean.
    y = y-mean(y);
    y = cumsum(y);         		% Integrate the above fluctuation time series ('y').
    for i = 1:size(DFA_x,1);				% 'DFA_x(i)' is the window size, which increases uniformly on a log10 scale!
        D = zeros(floor(size(y,1)/(DFA_x(i)*(1-DFA_Overlap))),1);		% initialize vector for temporarily storing the root-mean-square of each detrended window.
        tt = 0;
        for nn = 1:round(DFA_x(i)*(1-DFA_Overlap)):size(y,1)-DFA_x(i);	% we are going to detrend all windows in steps of 'n*(1-DFA_Overlap)'.
            tt=tt+1;
            D(tt) = (mean(fastdetrend(y(nn:nn+DFA_x(i))).^2,1))^(1/2);		% the square of the fluctuation around the local trend (of the window).
        end
        DFA_y(i) = mean(D(1:tt),1);						% the root-mean-square fluctuation of the integrated and detrended time series
    end  					  	       			% -- the F(n) in eq. (1) in Peng et al. 1995.
end

%% Supporting functions
function signal = fastdetrend(signal)
% A simple and fast detrend, see also the supporting function fastdetrend
% in the supporting functions folder
persistent a
persistent N
if (isempty(a) || size(signal,1) ~= N)
    N = size(signal,1);
    a = [zeros(N,1) ones(N,1)];
    a(1:N) = (1:N)'/N;
end
signal = signal - a*(a\signal);
end