function [xdat,ydat,ks] = power_plot(powerData,fig,noBins,marker,color)

    %this is to logbin 

    binn = 0:(max(powerData));
    [ys_h, xs_h] = hist(powerData, binn);
    [xss, yss] = Logbinning_v2(xs_h, ys_h, noBins);
    hold on;
    xdat = log10(xss);
    ydat = log10(yss/sum(yss));
    ok = ~(isnan(xdat) | isnan(ydat) |isinf(xdat) | isinf(ydat));
    plot(xdat(find(ok)),ydat(find(ok)),'linestyle', 'none','Marker', marker,'Color',color,'MarkerSize',10);
    xdat = xdat(find(ok));
    ydat = ydat(find(ok));

    minX = 1;
    maxX = 2500;

    avs = powerData;
    avs = avs(avs >= minX & avs <= maxX);
    t = 10.^(linspace(log10(minX),log10(maxX),noBins));
    for i = 1:noBins
        fb(i) = (1 - (minX/maxX)^0.5)^-1 * (1 - (minX/t(i))^0.5);
    end
    for i = 1:noBins
        ys(i) = nnz(avs<t(i));
    end
    ys = ys / length(avs);

    ks= 1+ sum(fb - ys)/noBins;

end
