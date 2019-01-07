%状态完全可测时，基于Lyapunov稳定性理论的MRAC
clear all; close all;

h=0.01; L=40/h; %数值积分步长和仿真步数
Ap=[0 1;-6 -7]; Bp=[0;8]; %对象参数
Am=[0 1;-10 -5]; Bm=[0;2]; %参考模型参数
Sz=size(Bp); n=Sz(1); m=Sz(2); %状态向量、输入维数

P=[3 1;1 1]; %正定矩阵
R1=1*eye(m); R2=1*eye(m); %自适应律参数矩阵

F0=zeros(m,n); K0=zeros(m); %初值
yr0=zeros(m,1); u0=zeros(m,1); e0=zeros(n,1);
xp0=zeros(n,1); xm0=zeros(n,1);

for k=1:L
    time(k)=k*h;
    yr(k)=1*sin(0.01*pi*time(k))+4*sin(0.2*pi*time(k))+sin(1*pi*time(k)); %输入信号
    xp(:,k)=xp0+h*(Ap*xp0+Bp*u0); %计算xp
    xm(:,k)=xm0+h*(Am*xm0+Bm*yr0); %计算xm
    e(:,k)=xm(:,k)-xp(:,k); %e=xm-xp
    
    F=F0+h*(R1*Bm'*P*e0*xp0'); %自适应律
    K=K0+h*(R2*Bm'*P*e0*yr0');
    u(:,k)=K*yr(k)+F*xp(:,k); %控制量
    
    parae(1,k)=norm(Am-Ap-Bp*F); %参数收敛于Am的偏差
    parae(2,k)=norm(Bm-Bp*K); %参数收敛于Bm的偏差
    
    %更新数据
    yr0=yr(:,k); u0=u(:,k); e0=e(:,k);
    xp0=xp(:,k); xm0=xm(:,k);
    F0=F; K0=K;
end
figure(1)
subplot(2,1,1);
plot(time,xm(1,:),'r',time,xp(1,:),':');
xlabel('t'); ylabel('x_m_1(t)、x_p_1(t)');
legend('x_m_1(t)','x_p_1(t)');
subplot(2,1,2);
plot(time,xm(2,:),'r',time,xp(2,:),':');
xlabel('t'); ylabel('x_m_2(t)、x_p_2(t)');
legend('x_m_2(t)','x_p_2(t)');
figure(2)
plot(time,parae(1,:),':',time,parae(2,:),'r');
xlabel('t'); ylabel('参数收敛偏差E');
legend('||A_m-A_p-B_p*F||_2','||B_m-B_p*K||_2');