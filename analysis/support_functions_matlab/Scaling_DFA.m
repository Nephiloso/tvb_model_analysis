function [DFA_exp,x,y] = Scaling_DFA(Data,Fs,DFA_SmallTime,DFA_LargeTime,DFA_SmallTimeFit,DFA_LargeTimeFit,DFA_Overlap,DFA_Plot,color,offset_y);
%
% Modified klaus.linkenkaer@cncr.vu.nl, 060502.
% Modified poil@get2net.dk, 17-april-2007 (added cumsum analysis)
% Modified poil@get2net.dk, 07-december-2007 (added improved plotting)
% Modified by S.-S. Poil, 07-feb-2008 (added improved detrending)
%******************************************************************************************************************
% Method references...
%
% Peng et al., Mosaic organization of DNA nucleotides, Phys. rev. E (49), 1685-1688 (1994).
% or for a better description: Peng et al., Quantification of scaling exponents and crossover phenomena
% in nonstationary heartbeat time series, Chaos (5), 82-87 (1995).
%
%
%******************************************************************************************************************
% Input parameters...
%
% Data	  		    : input data (vector).
% Fs			    : sampling frequency.
% DFA_SmallTime 	: minimum time-window size computed (in units of seconds!).
% DFA_LargeTime 	: maximum time-window size computed (in units of seconds!).
% DFA_SmallTimeFit	: smallest time scale (window size) to include in power-law fit (in units of seconds!).
% DFA_LargeTimeFit	: largest time scale (window size) to include in power-law fit (in units of seconds!).
% DFA_Overlap		: amount of DFA_Overlap between windows (to obtain higher SNR for the fluctuation estimates).
% DFA_Plot		    : should either be an axeshandle or an integer > 0
%
%
%******************************************************************************************************************
% output parameters...
%
% DFA_x			: Units in seconds. The logarithmically binned x-axis of window sizes (equidistant on a log axis).
% DFA_y			: the DFA fluctuation function, the y-axis.
% DFA_exp		: the DFA power-law exponent.
%
%
%******************************************************************************************************************
% Default parameters...

if size(Data,1) > 2
    Data=Data';
end

res_logbin = 10;	% number of bins pr decade, i.e., the spacing of the logarithmic scale.


%******************************************************************************************************************
% Defining window sizes to be log-linearly increasing.

d1 = floor(log10(DFA_SmallTime*Fs));
d2 = ceil(log10(DFA_LargeTime*Fs));
DFA_x_t = round(logspace(d1,d2,(d2-d1)*res_logbin));	% creates vector from 10^d1 to 10^d2 with N log-equidistant points.
DFA_x = DFA_x_t(find(DFA_SmallTime*Fs <= DFA_x_t & DFA_x_t <= DFA_LargeTime*Fs));	% Only include log-bins in the time range of interest!


%******************************************************************************************************************
% The DFA algorithm...

DFA_y = zeros(1,size(DFA_x,2));
y = zeros(1,size(Data,2));

b = Data./mean(Data);
b = b - mean(b);
% b = Data-mean(Data);		% First we convert the time series to a series of fluctuations b(i) around the mean.
%disp('NON_CUM_SUM');
y = b;
y = cumsum(b);         		% Integrate the above fluctuation time series ('b').
%% experimental msum stuff ,begin
%msum_i = 2:size(y,2);
%msum = mean(abs((y(msum_i)./y(msum_i-1)))); %cumsum analysis
%disp(msum)
%% experimental msum stuff ,end
for i = 1:size(DFA_x,2);				% 'DFA_x(i)' is the window size, which increases uniformly on a log10 scale!
    D = zeros(1,floor(size(y,2)/(DFA_x(i)*(1-DFA_Overlap))));		% initialize vector for temporarily storing the root-mean-square of each detrended window.
    tt = 0;
    for nn = 1:round(DFA_x(i)*(1-DFA_Overlap)):size(y,2)-DFA_x(i);	% we are going to detrend all windows in steps of 'n*(1-DFA_Overlap)'.
        tt=tt+1;
        D(tt) = (mean(fastdetrend(y(nn:nn+DFA_x(i))).^2))^(1/2);		% the square of the fluctuation around the local trend (of the window).
        %D(tt) = var(y(nn:nn+DFA_x(i)))^0.5;%        (mean((y(nn:nn+DFA_x(i))).^2))^(1/2);
       % disp('ERROR');

    end
    DFA_y(i) = mean(D(1:tt),2);						% the root-mean-square fluctuation of the integrated and detrended time series
end  					  	       			% -- the F(n) in eq. (1) in Peng et al. 1995.


%******************************************************************************************************************
% Fitting power-law and plotting...

DFA_SmallTime_LogSample = min(find(DFA_x>=DFA_SmallTime*Fs));		%
DFA_LargeTime_LogSample = max(find(DFA_x<=DFA_LargeTime*Fs));
DFA_SmallTimeFit_LogSample = min(find(DFA_x>=DFA_SmallTimeFit*Fs));
DFA_LargeTimeFit_LogSample = max(find(DFA_x<=DFA_LargeTimeFit*Fs));
X = [ones(1,DFA_LargeTimeFit_LogSample-DFA_SmallTimeFit_LogSample+1)' log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample))'];
Y = log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample))';

[DFA_exp,bint,r,rint,stats] = regress(Y,X);
DFA_exp = DFA_exp(2,1);
conf = ((bint(2,2))-(bint(2,1)))/2;		% compute +- 95% confidence intervals
%DFA_exp = X\Y; % for faster fitting ..
%DFA_exp  = DFA_exp(2);

if(DFA_Plot ~=0)
    %plot output
%     if ~ishandle(DFA_Plot)		%see if any figure handle is set
%        figure(DFA_Plot)
%        DFA_Plot = axes;
%     end

    %subplot(2,2,4)
%     axes(DFA_Plot)
    hold on
    
    x_all = log10(DFA_x/Fs);
    y_all = log10(DFA_y);
    
    
    x_fit = log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample)/Fs);
    y_fit = log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample));
    
    B = [ones(size(x_fit')), x_fit']\y_fit';
    Slope = B(2);
    disp(Slope);
    Intercept = B(1);
    
    x1 = x_fit(1);
    x2 = x_fit(end);
    y1 = x1*Slope+Intercept;
    y2 = x2*Slope+Intercept;
    
    y2 = y2-y1 + offset_y;
    y_all = y_all - y1 + offset_y;
    y1 = y1-y1 + offset_y;
    
    plot(x_all,y_all,'.','Color',color);
    
%     disp([x1 x2]);
    
    plot([x1 x2],[y1 y2],'Color',color,'LineWidth',2);
%     ll = refline([Slope Intercept]);
%     set(ll(1),'LineWidth',2);
    
  %  lsline_sisp
%     plot(log10(DFA_x(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)/Fs),log10(DFA_y(DFA_SmallTime_LogSample:DFA_LargeTime_LogSample)),strcat(color,'.'))
%     grid on
%     zoom on

    xlim([0 2]);
    ylim([0 1.6]);
    
    xlabel('log_{10}(time), [Seconds]','Fontsize',12)
    ylabel('log_{10} F(time)','Fontsize',12)
    
%     title(['DFA-exp=', num2str(DFA_exp), ' +-=', num2str(conf,3),', R^2 = ', num2str(stats(1,1),3)],'Fontsize',12)

end

%******************************************************************************************************************
% Conversions...

DFA_x = DFA_x/Fs;
end

function signal = fastdetrend(signal)
% fast detrending of "signal"
n = size(signal,1);
if n == 1,
    signal = signal(:); % make signal a row vector
end

% set up linear fitting 
N = size(signal,1);
a = [zeros(N,1) ones(N,1)];
a(1:N) = (1:N)'/N;

signal = signal - a*(a\signal); % remove best fit

if(n==1)
    signal = signal.'; % return correct dimensions
end
end

