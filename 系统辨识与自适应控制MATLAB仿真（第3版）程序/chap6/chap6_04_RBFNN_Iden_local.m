%基于局部误差的RBF神经网络辨识
clear all; close all;

ny=2; nu=3; d=2; %ny、nu、d为系统结构参数

L=600; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值

%设置RBF网络参数
n=ny+nu-d+1; m=6; %n、m分别为输入层和隐含层节点数
eta=0.5; %学习速率
alpha=0.05; %动量因子
ck1=20*ones(m,n); ck2=ck1; %隐含层中心向量初值: cki表示c(k-i)
bk1=40*ones(m,1); bk2=bk1; %隐含层中心向量初值: bki表示b(k-i)
wk1=rands(1,m); wk2=wk1; %隐含层至输出层权值的初值: wki表示w(k-i)
R=zeros(m,1); %定义R结构
db=zeros(m,1); %定义Δb结构
dc=zeros(m,n); %定义Δc结构

for k=1:L
    time(k)=k;
    u(k)=0.8*sin(0.01*pi*k); %控制输入信号
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
    %计算RBF网络输出
    x=[yk; uk(d:nu)];
    for j=1:m
        R(j)=exp(-norm(x-ck1(j,:)')^2/(2*bk1(j)^2));
    end 
    ym(k)=wk1*R;
    
    e(k)=y(k)-ym(k); %模型误差
    
    %RBF网络训练
    dw=eta*e(k)*R'; %Δw(k)
    for j=1:m
        db(j)=eta*e(k)*wk1(j)*R(j)*norm(x-ck1(j,:)')^2/bk1(j)^3; %Δb(k)
        for i=1:n
            dc(j,i)=eta*e(k)*wk1(j)*R(j)*(x(i)-ck1(j,i))/bk1(j)^2; %Δc(k)
        end
    end

    w=wk1+dw+alpha*(wk1-wk2);
    b=bk1+db+alpha*(bk1-bk2);
    c=ck1+dc+alpha*(ck1-ck2);
    
    %更新数据
    bk2=bk1; bk1=b;
    ck2=ck1; ck1=c;
    wk2=wk1; wk1=w;
    
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
legend('y(k)','y_m(k)'); axis([0 L -.5 2]);
subplot(212)
plot(time,y-ym,'k');
xlabel('k'); ylabel('e(k)'); %axis([0 L -.5 .5]);