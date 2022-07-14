load('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\grid_search\results.mat')
load('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\grid_search\results2.mat')
T0=T;
load('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\grid_search\results.mat')
T_tot = [T0;T];

c_ee = 11:1:16;
c_ei = 2:1:15;
[X,Y] = meshgrid(c_ee,c_ei);
Z = zeros(length(c_ei),length(c_ee));
for i=c_ee
    for j=c_ei
        idx = find(T_tot.c_ee==i&T_tot.c_ei==j);
        if idx
            Z(j-1,i-10) = T_tot.DFA_exp_delta(idx);
        end
        clear idx
    end
end
contour(X,Y,Z)
xlabel('W[E->E]')
ylabel('W[I->E]')
colorbar