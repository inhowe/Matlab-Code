%基于局部误差的BP神经网络辨识
clear all; close all;

ny=2; nu=3; d=2; %ny、nu、d为系统结构参数

L=600; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值

%设置BP网络参数
n=ny+nu-d+1; m=6; %n、m分别为输入层和隐含层节点数
eta=0.5; %学习速率
alpha=0.05; %动量因子
w1k1=[-0.8633    0.9231    0.1046    0.7128
   -0.1273    0.5248   -0.5638   -0.1951
   -0.6523   -0.9853    0.5447   -0.3640
   -0.9478    0.3601   -0.5439    0.2173
    0.9094    0.4119   -0.2583    0.8204
   -0.1388    0.2903    0.7819    0.8182];
%w1k1=rands(m,n); %输入层至隐含层权值的初值: w1ki表示w1(k-i)
w1k2=w1k1;
w2k1=[0.1832  -0.3349  0.7061  -0.1152  0.8087  -0.9336];
%w2k1=rands(1,m); %隐含层至输出层权值的初值: w2ki表示w2(k-i)
w2k2=w2k1;

for  k=1:L
    time(k)=k;
    u(k)=0.8*sin(0.01*pi*k); %控制输入信号
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
    %计算BP网络输出
    X=[yk; uk(d:nu)];
    O1=X; 
    net2=w1k1*O1;
    O2=1./(1+exp(-net2));    
    ym(k)=w2k1*O2;
    
    e(k)=y(k)-ym(k); %模型误差

    %BP网络训练
    dw2=eta*e(k)*O2'; 
    w2=w2k1+dw2+alpha*(w2k1-w2k2); %w2(k)

    df=exp(-net2)./(1+exp(-net2)).^2; %激励函数的导数
    dw1=eta*e(k)*w2k1'.*df*O1'; %矩阵形式计算
    %for j=1:m 
        %for i=1:n
            %dw1(j,i)=eta*e(k)*w2k1(j)*df(j)*O1(i); %标量形式计算
        %end
    %end
    w1=w1k1+dw1+alpha*(w1k1-w1k2); %w1(k)
    
    %更新数据
    w1k2=w1k1; w1k1=w1;
    w2k2=w2k1; w2k1=w2;
    
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
plot(time,y,'r:',time,ym,'k');
xlabel('k'); ylabel('y(k)、y_m(k)');
legend('y(k)','y_m(k)'); %axis([0 L -.4 1.6]);
subplot(212)
plot(time,y-ym,'k');
xlabel('k'); ylabel('e(k)'); axis([0 L -1 1]);