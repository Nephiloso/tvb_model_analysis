%The script aims to plot all the basic things as Klaus's paper in 2004.
dt=1;
Fs=1000/dt;
Fs = %1000;  % monitors-period: ms => Fs: Hz
result = readmatrix('C:\Users\wpp_1\Downloads\2022-04-19_simu1_results.csv');  % ! heading problem
time = result(:,1)/1000;
v1 = result(:,2);

plot_time_intv=[100,800]; % ms
plot_intv= [double2int(plot_time_intv(1)/dt):double2int(plot_time_intv(2)/dt)];
figure();
plot(time[plot_intv],v1[plot_inv])