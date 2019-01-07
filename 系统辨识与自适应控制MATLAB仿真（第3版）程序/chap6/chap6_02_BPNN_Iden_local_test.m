%基于局部误差的BP神经网络辨识（模型测试）
clear all; close all;

ny=2; nu=3; d=2; %ny、nu、d为系统结构参数

L=600; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值

%BP网络训练后获得的参数（k=600）
w1=[-0.8636    0.9063    0.1142    0.7152
   -0.1340    0.5349   -0.5624   -0.1881
   -0.6303   -1.0082    0.5633   -0.3626
   -0.9515    0.3592   -0.5399    0.2215
    0.9268    0.3838   -0.2346    0.8257
   -0.1506    0.3201    0.7703    0.8224]; %输入层至隐含层权值
w2=[0.3360  -0.1615  0.9604  0.0587  1.0013  -0.7268]; %隐含层至输出层权值

%BP网络模型测试
for  k=1:2*L
    time(k)=k;
    u(k)=0.8*sin(0.01*pi*k); %控制输入信号
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
    %计算BP网络输出
    X=[yk; uk(d:nu)];
    O1=X;     
    net2=w1*O1;
    O2=1./(1+exp(-net2));   
    ymt(k)=w2*O2;
    
    et(k)=y(k)-ymt(k); %模型误差
    
    %更新数据
    for i=nu:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=ny:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
subplot(211)
plot(time,y,'r:',time,ymt,'k');
xlabel('k'); ylabel('y(k)、y_m(k)');
legend('y(k)','y_m(k)'); %axis([0 L -.4 1.6]);
subplot(212)
plot(time,et,'k');
xlabel('k'); ylabel('e(k)'); axis([0 L -1.5 1]);