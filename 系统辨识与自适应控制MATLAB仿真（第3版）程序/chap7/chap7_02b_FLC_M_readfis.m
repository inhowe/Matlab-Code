%非线性系统模糊逻辑控制（readfis + 实时控制）
clear all; close all;

a=readfis('fuzcon_m'); %从磁盘当前目录中装入已建好的FIS

ny=2; nu=3; d=2; Ts=1;  %ny、nu、d为系统结构参数,Ts为采样周期

L=1500; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值
ek1=0; %e(k-1)

ke=12; kec=12; kuc=0.005; %量化因子
for k=1:L
    time(k)=k*Ts;
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
    yr(k)=0.25*sign(sin(0.002*pi*k))+0.75; %采集期望输出数据
    e(k)=yr(k)-y(k);
    ec=(e(k)-ek1)/Ts;
    
    e1=ke*e(k); ec1=kec*ec; %模糊控制器输入
    if e1>6
        e1=6;
    end
    if e1<-6
        e1=-6;
    end
    if ec1>6
        ec1=6;
    end
    if ec1<-6
        ec1=-6;
    end

    uc=kuc*evalfis([e1 ec1],a); %执行模糊推理计算
    u(k)=uk(1)+uc;
    
    %更新数据
    for i=nu:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=ny:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
    
    ek1 = e(k); 
end
figure(1)
plot(time,yr,'r--',time,y,'k');
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L 0.4 1.1]);
figure(2);
plot(time,u,'k');
xlabel('k'); ylabel('u(k)'); axis([0 L -.8 .8]);