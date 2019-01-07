%可调增益MIT-MRAC归一化算法
clear all; close all;

h=0.1; L=100/h; %数值积分步长和仿真步数
num=[1]; den=[1 1 1]; n=length(den)-1; %对象参数
kp=1; [Ap,Bp,Cp,Dp]=tf2ss(kp*num,den); %传递函数型转换为状态空间型
km=1; [Am,Bm,Cm,Dm]=tf2ss(km*num,den); %参考模型参数

gamma=0.1; %自适应增益
alpha=0.01; beta=2;

yr0=0; u0=0; e0=0; ym0=0; %初值
xp0=zeros(n,1); xm0=zeros(n,1); %状态向量初值
kc0=0; %可调增益初值
r=1.2; yr=r*[ones(1,L/4) -ones(1,L/4) ones(1,L/4) -ones(1,L/4)]; %输入信号

for k=1:L
    time(k)=k*h;
    xp(:,k)=xp0+h*(Ap*xp0+Bp*u0); 
    yp(k)=Cp*xp(:,k)+Dp*u0; %计算yp
    
    xm(:,k)=xm0+h*(Am*xm0+Bm*yr0); 
    ym(k)=Cm*xm(:,k)+Dm*yr0; %计算ym
    e(k)=ym(k)-yp(k); %e=ym-yp
    
    DD=e0*ym0/km/(alpha+(ym0/km)^2);

    kc=kc0+h*gamma*DD; %MIT自适应律
    u(k)=kc*yr(k); %控制量
    
    %更新数据
    yr0=yr(k); u0=u(k); e0=e(k); ym0=ym(k);
    xp0=xp(:,k); xm0=xm(:,k);
    kc0=kc;
end
plot(time,ym,'r',time,yp,':');
xlabel('t'); ylabel('y_m(t)、y_p(t)');
legend('y_m(t)','y_p(t)');