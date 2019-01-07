%单步Diophantine方程的求解
clear all;

a=[1 -1.7 0.7]; b=[1 0.5]; c=[1 0.2]; d=4; %例4.1（1）
%a=[1 -1.7 0.7]; b=[1 2]; c=[1 0.2]; d=4; %例4.1（2）

[e,f,g]=sindiophantine(a,b,c,d) %调用函数sindiophantine（结果见MATLAB命令窗口）