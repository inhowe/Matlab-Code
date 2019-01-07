%增量式PID控制
clear all; close all;

den=[1 1 0]; num=[1]; Ts=0.5; Td=1; %连续系统对象
sys=tf(num,den,'inputdelay',Td);
dsys=c2d(sys,Ts,'zoh');
[dnum,a]=tfdata(dsys,'v');
na=length(a)-1; b=dnum(2:na+1); nb=length(b)-1;
d=Td/Ts+1;

kp=0.4; ki=0.0; kd=1; %PID控制器参数（试凑法）

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
ek=zeros(2,1); %输出误差初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4,1)]; %期望输出

for k=1:L
    time(k)=k*Ts;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb); %采集输出数据
    e(k)=yr(k)-y(k);
    
    %增量式PID控制律
    du=kp*(e(k)-ek(1))+ki*e(k)+kd*(e(k)-2*ek(1)+ek(2));
    u(k)=uk(1)+du;

    %更新数据
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
    ek(2)=ek(1);ek(1)=e(k);
end
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('t'); ylabel('y_r(t)、y(t)');
legend('y_r(t)','y(t)');
subplot(2,1,2);
plot(time,u);
xlabel('t'); ylabel('u(t)');