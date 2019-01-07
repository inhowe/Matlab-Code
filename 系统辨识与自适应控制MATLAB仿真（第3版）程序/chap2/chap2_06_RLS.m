%递推最小二乘参数估计（RLS）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3; %对象参数
na=length(a)-1; nb=length(b)-1; %na、nb为A、B阶次

L=400; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
u=randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

theta=[a(2:na+1);b]; %对象参数真值

thetae_1=zeros(na+nb+1,1); %thetae初值
P=10^6*eye(na+nb+1); 
for k=1:L
    % 被辨识的模型，是迭代生成的形式
    phi=[-yk;uk(d:d+nb)]; %此处phi为列向量
    y(k)=phi'*theta+xi(k); %采集输出数据
   
    %递推最小二乘法
    K=P*phi/(1+phi'*P*phi);
    thetae(:,k)=thetae_1+K*(y(k)-phi'*thetae_1);
    P=(eye(na+nb+1)-K*phi')*P;
    
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
plot([1:L],thetae); %line([1,L],[theta,theta]);
xlabel('k'); ylabel('参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 1.5]);