%递推极大似然参数估计（RML）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; c=[1 -0.5]'; d=1; %对象参数
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为A、B、C阶次
nn=max(na,nc); %用于yf(k-i)、uf(k-i)更新

L=1000; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
xik=zeros(nc,1); %白噪声初值
xiek=zeros(nc,1); %白噪声估计初值
yfk=zeros(nn,1); %yf(k-i)
ufk=zeros(nn,1); %uf(k-i)
xiefk=zeros(nc,1); %ξf(k-i)
u=randn(L,1); %输入采用白噪声序列
xi=randn(L,1); %白噪声序列

thetae_1=zeros(na+nb+1+nc,1); %参数估计初值
P=eye(na+nb+1+nc);
for k=1:L
    y(k)=-a(2:na+1)'*yk+b'*uk(d:d+nb)+c'*[xi(k);xik]; %采集输出数据
        
    %构造向量
    phi=[-yk;uk(d:d+nb);xiek];
    xie=y(k)-phi'*thetae_1;
    phif=[-yfk(1:na);ufk(d:d+nb);xiefk];
    
    %递推极大似然参数估计算法
    K=P*phif/(1+phif'*P*phif);
    thetae(:,k)=thetae_1+K*xie;
    P=(eye(na+nb+1+nc)-K*phif')*P;    
        
    yf=y(k)-thetae(na+nb+2:na+nb+1+nc,k)'*yfk(1:nc); %yf(k)
    uf=u(k)-thetae(na+nb+2:na+nb+1+nc,k)'*ufk(1:nc); %uf(k)
    xief=xie-thetae(na+nb+2:na+nb+1+nc,k)'*xiefk(1:nc); %xief(k)
      
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
        xiefk(i)=xiefk(i-1);
    end
    xik(1)=xi(k);
    xiek(1)=xie;
    xiefk(1)=xief;
    
    for i=nn:-1:2
        yfk(i)=yfk(i-1);
        ufk(i)=ufk(i-1);
    end
    yfk(1)=yf;
    ufk(1)=uf;
end
figure(1)
plot([1:L],thetae(1:na,:),[1:L],thetae(na+nb+2:na+nb+1+nc,:));
xlabel('k'); ylabel('参数估计a、c');
legend('a_1','a_2','c_1'); axis([0 L -2 2]);
figure(2)
plot([1:L],thetae(na+1:na+nb+1,:));
xlabel('k'); ylabel('参数估计b');
legend('b_0','b_1'); axis([0 L 0 1.5]);