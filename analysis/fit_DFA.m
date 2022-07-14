function DFA_exp = fitDFA( DFA_x, DFA_y, FitInterval,Fs )
    if isempty(DFA_x) || isempty(DFA_y)
        DFA_exp = nan;
        return
    end
    if isempty(Fs)
        Fs=1000;
    end
    DFA_SmallTimeFit_LogSample = min(find(DFA_x>=FitInterval(1)*Fs));
    DFA_LargeTimeFit_LogSample = max(find(DFA_x<=FitInterval(2)*Fs));
    X = [ones(1,DFA_LargeTimeFit_LogSample-DFA_SmallTimeFit_LogSample+1)' log10(DFA_x(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample))'];
    Y = log10(DFA_y(DFA_SmallTimeFit_LogSample:DFA_LargeTimeFit_LogSample));
    DFA_exp = X\Y; %least-square fit
    DFA_exp = DFA_exp(2);
end

