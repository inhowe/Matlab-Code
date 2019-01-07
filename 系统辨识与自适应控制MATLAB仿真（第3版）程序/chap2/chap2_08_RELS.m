%递推增广最小二乘参数估计（RELS）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; c=[1 -1 0.2]'; d=3; %对象参数
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为A、B、C阶次

L=1000; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
xik=zeros(nc,1); %噪声初值
xiek=zeros(nc,1); %噪声估计初值
u=randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

theta=[a(2:na+1);b;c(2:nc+1)]; %对象参数

thetae_1=zeros(na+nb+1+nc,1); %na+nb+1+nc为辨识参数个数
P=10^6*eye(na+nb+1+nc);
for k=1:L
    phi=[-yk;uk(d:d+nb);xik];
    y(k)=phi'*theta+xi(k); %采集输出数据
    
    phie=[-yk;uk(d:d+nb);xiek]; %组建phie
    
    %递推增广最小二乘法
    K=P*phie/(1+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1+nc)-K*phie')*P;
    
    xie=y(k)-phie'*thetae(:,k); %白噪声的估计值
    
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
    
    for i=nc:-1:2
        xik(i)=xik(i-1);
        xiek(i)=xiek(i-1);
    end
    xik(1)=xi(k);
    xiek(1)=xie;
end
figure(1)
plot([1:L],thetae(1:na,:));
xlabel('k'); ylabel('参数估计a');
legend('a_1','a_2'); axis([0 L -2 2]);
figure(2)
plot([1:L],thetae(na+1:na+nb+1,:));
xlabel('k'); ylabel('参数估计b');
legend('b_0','b_1'); axis([0 L 0 1.5]);
figure(3)
plot([1:L],thetae(na+nb+2:na+nb+nc+1,:));
xlabel('k'); ylabel('参数估计c');
legend('c_1','c_2'); axis([0 L -2 2]);