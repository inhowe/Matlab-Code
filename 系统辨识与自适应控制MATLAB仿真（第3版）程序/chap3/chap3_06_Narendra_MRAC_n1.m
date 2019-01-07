%n-m=1 Narendra自适应控制律（不稳定对象）
clear all; close all;

h=0.01; L=40/h; %数值积分步长、仿真步数
nump=[1 1]; denp=[1 -5 6]; [Ap,Bp,Cp,Dp]=tf2ss(nump,denp); n=length(denp)-1; %对象参数
numm=[1 2]; denm=[1 3 6]; [Am,Bm,Cm,Dm]=tf2ss(numm,denm); %参考模型参数

Df=numm; %辅助信号发生器传递函数分母多项式
Af=[[zeros(n-2,1),eye(n-2)];-Df(n:-1:2)]; %辅助信号发生器状态矩阵
Bf=[zeros(n-2,1);1]; %辅助信号发生器输入矩阵

yr0=0; yp0=0; u0=0; e0=0; %初值
v10=zeros(n-1,1); v20=zeros(n-1,1); %辅助信号发生器状态初值
xp0=zeros(n,1); xm0=zeros(n,1); %状态向量初值
theta0=zeros(2*n,1); %可调参数向量初值
r=2; yr=r*[ones(1,L/4) -ones(1,L/4) ones(1,L/4) -ones(1,L/4)]; %参考输入信号

Gamma=10*eye(2*n); %自适应增益矩阵（正定矩阵）
for k=1:L
    time(k)=k*h;
    xp(:,k)=xp0+h*(Ap*xp0+Bp*u0);
    yp(k)=Cp*xp(:,k)+Dp*u0; %计算yp
    
    xm(:,k)=xm0+h*(Am*xm0+Bm*yr0);
    ym(k)=Cm*xm(:,k)+Dm*yr0; %计算ym
    e(k)=ym(k)-yp(k); %e=ym-yp
    
    v1=v10+h*(Af*v10+Bf*u0); %计算v1
    v2=v20+h*(Af*v20+Bf*yp0); %计算v1
    
    phi0=[yr0;v10;yp0;v20]; %组建k-1时刻的数据向量
    theta(:,k)=theta0+h*e0*Gamma*phi0; %自适应律
    phi=[yr(k);v1;yp(k);v2]; %组建k时刻的数据向量
    u(k)=theta(:,k)'*phi; %自适应控律
    
    %更新数据
    yr0=yr(k); yp0=yp(k); u0=u(k); e0=e(k);
    v10=v1; v20=v2;
    xp0=xp(:,k); xm0=xm(:,k);
    phi0=phi; theta0=theta(:,k);
end
subplot(2,1,1);
plot(time,ym,'r',time,yp,':');
xlabel('t'); ylabel('y_m(t)、y_p(t)');
legend('y_m(t)','y_p(t)');
subplot(2,1,2);
plot(time,u);
xlabel('t'); ylabel('u(t)');