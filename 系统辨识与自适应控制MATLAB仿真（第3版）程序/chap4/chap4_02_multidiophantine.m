%多步Diophantine方程的求解
clear all;

a=[1  -3  3.1  -1.1]; b=[1 2]; c=[1];
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %A、B、C的阶次
N=2; %预测步数

[E,F,G]=multidiophantine(a,b,c,N) %调用函数multidiophantine（结果见MATLAB命令窗口）