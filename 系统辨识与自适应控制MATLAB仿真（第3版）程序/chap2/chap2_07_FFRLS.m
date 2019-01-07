%遗忘因子递推最小二乘参数估计（FFRLS）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3; %对象参数
na=length(a)-1; nb=length(b)-1; %na、nb为A、B阶次

L=1000; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
u=randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

thetae_1=zeros(na+nb+1,1); %thetae初值
P=10^6*eye(na+nb+1);
lambda=1; %遗忘因子范围[0.9 1]
for k=1:L
    if k==501
        a=[1 -1 0.4]';b=[1.5 0.2]'; %对象参数突变
    end
    theta(:,k)=[a(2:na+1);b]; %对象参数真值
    
    phi=[-yk;uk(d:d+nb)];
    y(k)=phi'*theta(:,k)+xi(k); %采集输出数据
   
    %遗忘因子递推最小二乘法
    K=P*phi/(lambda+phi'*P*phi);
    thetae(:,k)=thetae_1+K*(y(k)-phi'*thetae_1);
    P=(eye(na+nb+1)-K*phi')*P/lambda;
    
    %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
subplot(1,2,1)
plot([1:L],thetae(1:na,:)); hold on; plot([1:L],theta(1:na,:),'k:');
xlabel('k'); ylabel('参数估计a');
legend('a_1','a_2'); axis([0 L -2 2]);
subplot(1,2,2)
plot([1:L],thetae(na+1:na+nb+1,:)); hold on; plot([1:L],theta(na+1:na+nb+1,:),'k:');
xlabel('k'); ylabel('参数估计b');
legend('b_0','b_1'); axis([0 L -0.5 2]);