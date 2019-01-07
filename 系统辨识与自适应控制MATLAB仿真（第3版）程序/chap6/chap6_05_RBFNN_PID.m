%基于局部误差的RBF神经网络辨识+PID自校正控制
clear all; close all;

ny=2; nu=3; d=2; %ny、nu、d为系统结构参数

L=1500; %仿真长度
uk=zeros(nu,1); %控制输入初值：uk(i)表示u(k-i);
yk=zeros(ny,1); %系统输出初值

%设置RBF网络参数
n=ny+nu+1; m=10; %n、m分别为输入层和隐含层节点数（注意：输入层和隐含层节点数与辨识程序chap6_04不同！）
eta=0.5; %学习速率
alpha=0.05; %动量因子
ck1=20*ones(m,n); ck2=ck1; %隐含层中心向量初值: cki表示c(k-i)
bk1=40*ones(m,1); bk2=bk1; %隐含层中心向量初值: bki表示b(k-i)
wk1=[0.3219  0.4595  0.7815  0.9646  0.5381  0.1629  0.8566  0.1602  -0.9660  -0.7583];
%wk1=rands(1,m);
wk2=wk1; %隐含层至输出层权值的初值: wki表示w(k-i)
R=zeros(m,1); %定义R结构
db=zeros(m,1); %定义Δb结构
dc=zeros(m,n); %定义Δc结构

%设置PID参数初值
Kp1=0.0; Kp2=Kp1; %比例: Kpi表示Kp(k-i)
Ki1=0.0; Ki2=Ki1; %积分: Kii表示Ki(k-i)
Kd1=0.0; Kd2=Kd1; %微分: Kdi表示Kd(k-i)
eck1=0;  eck2=0; %误差: ecki表示ec(k-i)

etac=1; %PID参数学习速率
alphac=0.1; %PID参数动量因子

for k=1:L
    time(k)=k;
    y(k)=uk(2)^3+uk(3)^3+(0.8+yk(1)^3)/(1+yk(1)^2+yk(2)^4); %采集系统输出数据
    
    %计算PID控制量u(k)
    yr(k)=0.25*sign(sin(0.002*pi*k))+0.75;
    ec(k)=yr(k)-y(k);
    xc(1)=ec(k)-eck1;
    xc(2)=ec(k);
    xc(3)=ec(k)-2*eck1+eck2;
    
    du=Kp1*xc(1)+Ki1*xc(2)+Kd1*xc(3); %控制增量Δu(k)
    u(k)=uk(1)+du;
    
    %计算RBF网络输出
    x=[yk; u(k); uk]; %RBF网络输入（包含u(k)！！！） 
    for j=1:m
        R(j)=exp(-norm(x-ck1(j,:)')^2/(2*bk1(j)^2));
    end 
    ym(k)=wk1*R; 
    
    %计算Jacobian信息
    J(k)=0;
    for j=1:m
        J(k)=J(k)+wk1(j)*R(j)*(ck1(j,ny+1)-u(k))/bk1(j)^2;
    end
    
    %PID参数学习
    dKp=etac*ec(k)*J(k)*xc(1); %ΔKp(k)
    dKi=etac*ec(k)*J(k)*xc(2); 
    dKd=etac*ec(k)*J(k)*xc(3); 
    
    Kp(k)=Kp1+dKp+alphac*(Kp1-Kp2); %Kp(k)
    Ki(k)=Ki1+dKi+alphac*(Ki1-Ki2);
    Kd(k)=Kd1+dKd+alphac*(Kd1-Kd2);
    if Kp(k)<0
        Kp(k)=0;
    end
    if Ki(k)<0
        Ki(k)=0;
    end
    if Kd(k)<0
        Kd(k)=0;
    end
    
    %RBF网络训练
    e(k)=y(k)-ym(k); %模型误差
    
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
    
    eck2=eck1; eck1=ec(k);
    Kp2=Kp1; Kp1=Kp(k);
    Ki2=Ki1; Ki1=Ki(k);
    Kd2=Kd1; Kd1=Kd(k);
end
figure(1);
plot(time,yr,'r--',time,y,'k:',time,ym,'k');
xlabel('k'); ylabel('y_r(k)、y(k)、y_m(k)');
legend('y_r(k)','y(k)','y_m(k)'); axis([0 L 0.4 1.1]);
figure(2);
plot(time,y-ym,'k');
xlabel('k'); ylabel('e(k)'); axis([0 L -.1 .1]);
figure(3);
plot(time,u,'k');
xlabel('k'); ylabel('u(k)'); axis([0 L -.8 .8]);
figure(4);
subplot(311)
plot(time,Kp,'k');
set(gca,'xtick',[]); ylabel('Kp(k)');
subplot(312)
plot(time,Ki,'k');
set(gca,'xtick',[]); ylabel('Ki(k)'); 
subplot(313)
plot(time,Kd,'k');
xlabel('k'); ylabel('Kd(k)'); axis([0 L -.5 .5]);
figure(5)
plot(time,J,'k');
xlabel('k'); ylabel('dy(k)/du(k)'); 