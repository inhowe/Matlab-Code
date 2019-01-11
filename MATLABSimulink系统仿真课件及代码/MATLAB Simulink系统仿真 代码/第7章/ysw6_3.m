%PID 控制器
clc % 清屏
clear all; % 删除workplace变量
close all; % 关掉显示图形窗口
M=0.5;m=0.5;b=0.1;I=0.006;l=0.3;g=9.8;
a=(M+m)*m*g*l/((M+m)*I+M*m*l^2);b=-m*l/(((M+m)*I+M*m*l^2));
c=-m^2*l^2*g/((M+m)*I+M*m*l^2);d=(I+m*l^2)/((M+m)*I+M*m*l^2);
A=[           0                   1 0             0;...
    (M+m)*m*g*l/((M+m)*I+M*m*l^2) 0 0 m*l*b/((M+m)*I+M*m*l^2);...
              0                   0 0             1;...
     -m^2*l^2*g/((M+m)*I+M*m*l^2) 0 0 -(I+m*l^2)*b/((M+m)*I+M*m*l^2)];
B=[0;-m*l/(((M+m)*I+M*m*l^2));0;(I+m*l^2)/((M+m)*I+M*m*l^2)];
C=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
D=[0;0;0;0];
p2=eig(A)';             % A特征值求解
p=[-10,-7,-1.901,-1.9]; % 极点配置
K=place(A,B,p)          % 状态反馈矩阵
eig(A-B*K)'  % 极点逆向求解
%%仿真结果验证
[x,y]=sim('pedulumpid.mdl');
subplot(121),plot(y(:,1),'r','linewidth',2);
grid on,title('倾角控制')
subplot(122),plot(y(:,3),'r','linewidth',2);
grid on,title('位移控制')
