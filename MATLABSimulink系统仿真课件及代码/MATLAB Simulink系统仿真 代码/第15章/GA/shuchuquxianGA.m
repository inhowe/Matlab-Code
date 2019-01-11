%基于GA算法的PID参数优化
clc % 清屏
clear all; % 删除workplace变量
close all; % 关掉显示图形窗口
ts=0.001;
sys=tf([1.6],[1 1.5 1.6],'inputdelay',0.1);
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');
sim('ysw_PID1.slx');
figure(1)
time = 0:1/(length(simout)-1):1;
plot(time,1-simout,'b','LineWidth',2)
xlabel('time(s)'),ylabel('yout');
grid on
title('GA优化阶跃响应输出曲线')

figure(2)
plot(time,simout,'r--','LineWidth',2)
xlabel('time(s)'),ylabel('error');
grid on
title('GA优化阶跃响应输出误差曲线')


 