clear all
M= readtable('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures\critical_grid_analysis.csv');
delta=M.dfa_all_raw;
c_ee = 6:4:23;
c_ei = 6:4:23;

[X,Y] = meshgrid(c_ei,c_ee);
Z = zeros(length(c_ei),length(c_ee));
for i=c_ee
    for j=c_ei
        idx = find(M.config_c_ee==i&M.config_c_ei==j);
        if idx
            Z(c_ee==i,c_ei==j) = delta(idx);
        end
        clear idx
    end
end
%contour(X,Y,Z)
figure('color','w');
yvalues = cellstr (string(uint8(c_ee)));
xvalues = cellstr (string(uint8(c_ei)));
hm = heatmap(xvalues,yvalues,Z);
xlabel('W[I->E]')
ylabel('W[E->E]')
title('DFA-raw signal');
hm.YDisplayData = flipud(hm.YDisplayData);
colormap('jet');
caxis([0.5 1]);
colorbar