%基于全局误差的BP神经网络辨识
clear all; close all;

ny=2; nu=3; d=2; %ny、nu、d为系统结构参数

L=600; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值

%设置BP网络参数
n=ny+nu-d+1; m=6; %n、m分别为输入层和隐含层节点数
eta=0.5; %学习速率
alpha=0.05; %动量因子
%w1k1=rands(m,n); %输入层至隐含层权值的初值: w1ki表示w1(k-i)
w1k1=[ -0.8633    0.9231    0.1046    0.7128
   -0.1273    0.5248   -0.5638   -0.1951
   -0.6523   -0.9853    0.5447   -0.3640
   -0.9478    0.3601   -0.5439    0.2173
    0.9094    0.4119   -0.2583    0.8204
   -0.1388    0.2903    0.7819    0.8182];
w1k2=w1k1; 
%w2k1=rands(1,m); %隐含层至输出层权值的初值: w2ki表示w2(k-i)
w2k1=[0.1832  -0.3349  0.7061  -0.1152  0.8087  -0.9336];
w2k2=w2k1; 

%产生训练样本2.5L个
for  k=1:2.5*L
    time(k)=k;
    u(k)=0.8*sin(0.01*pi*k); %控制输入信号
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
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

%利用X(k)训练BP网络，k=100, 101, ..., 900
eg=10; %初始化全局误差
epsilon=0.002; %全局误差精度
num=0; %初始化训练步数
M=100; %最大训练次数

while(eg>epsilon) %由全局误差控制训练次数
%while(num<M) %直接设定训练次数
    num=num+1;
    es(num)=0;
    for k=100:1.5*L        
        %计算BP网络输出
        X=[y(k-1); y(k-2); u(k-2); u(k-3)];
        O1=X;    
        net2=w1k1*O1;
        O2=1./(1+exp(-net2));    
        ym(k)=w2k1*O2; 
    
        e(k)=y(k)-ym(k); %模型误差
        es(num)=es(num)+e(k)^2/2; %累计误差平方   

        %BP网络训练
        dw2=eta*e(k)*O2'; 
        w2=w2k1+dw2+alpha*(w2k1-w2k2); %w2(k)

        df=exp(-net2)./(1+exp(-net2)).^2; %激励函数的导数
        dw1=eta*e(k)*w2k1'.*df*O1'; %矩阵形式计算
        w1=w1k1+dw1+alpha*(w1k1-w1k2); %w1(k)
    
        %更新数据
        w1k2=w1k1; w1k1=w1;
        w2k2=w2k1; w2k1=w2;
    end
    eg=es(num);
end %获得满足全局误差要求的权值w1、w2

%利用X(k)测试BP网络模型（即w1、w2），k=901, 902, ..., 1500
egt=0; %测试样本的输出全局误差
for k=1.5*L+1:2.5*L
    %计算BP网络输出
    X=[y(k-1); y(k-2); u(k-2); u(k-3)];
    O1=X;    
    net2=w1k1*O1;
    O2=1./(1+exp(-net2));   
    ymt(k)=w2k1*O2;
    
    et(k)=y(k)-ymt(k); %模型误差
    egt=egt+et(k)^2/2;
end

t1=100:1.5*L;
figure(1);
subplot(211)
plot(t1,y(t1),'r:',t1,ym(t1),'k');
xlabel('k'); ylabel('y(k)、y_m(k)');
legend('y(k)','y_m(k)'); %axis([0 L -.4 1.6]);
subplot(212)
plot(1:num,es,'k');
xlabel('Steps'); ylabel('E=\Sigma{e^2(k)/2}'); axis([0 num 0 max(es)]);
if num>500
    axes('Position',[0.3,0.25,0.4,0.16]); %生成子图
    t0=1:100;
    plot(t0,es(t0),'b'); axis([0 max(t0) 0 max(es)]);
end

t2=1.5*L+1:2.5*L;
figure(2);
subplot(211)
plot(t2,y(t2),'r:',t2,ymt(t2),'k');
xlabel('k'); ylabel('y(k)、y_m(k)');
legend('y(k)','y_m(k)'); %axis([0 L -.4 1.6]);
subplot(212)
plot(t2,et(t2),'k');
xlabel('k'); ylabel('e(k)'); %axis([0 L -.5 .5]);