%数值积分方法求解连续系统输出（龙格-库塔法与欧拉法的比较）
clear all; close all;

h=0.1; L=15/h; %h计算步长、L为仿真步数
z=[-1 -2]; p=[-4 -0.5+j -0.5-j]; k=2.5; %对象的零极点型
[A,B,C,D]=zp2ss(z,p,k); %转化为状态空间型
u=1*ones(L,1); u0=0; %输入及初值
n=length(p); %对象阶次

x0=zeros(n,1); xe0=zeros(n,1); %x0、xe0分别为龙格-库塔法、欧拉法的状态初值
for i=1:L
    time(i)=i*h;

    %欧拉法
    xe=xe0+h*(A*xe0+B*u0);
    ye(i)=C*xe;
    
    %龙格-库塔法
    k1=A*x0+B*u0;
    k2=A*(x0+h*k1/2)+B*u0;
    k3=A*(x0+h*k2/2)+B*u0;
    k4=A*(x0+h*k3)+B*u(i);
    x=x0+h*(k1+2*k2+2*k3+k4)/6;
    y(i)=C*x;
    
    %更新数据
    u0=u(i);x0=x;xe0=xe;
end
plot(time,u,'k-.',time,ye,':',time,y,'r');
xlabel('t'); ylabel('y_r(t)、y(t)');
legend('y_r(t)','Euler:y(t)','Runge-Kutta:y(t)');