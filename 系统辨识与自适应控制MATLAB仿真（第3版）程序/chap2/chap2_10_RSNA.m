%递推随机牛顿参数估计（RSNA）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3; %对象参数
na=length(a)-1; nb=length(b)-1; %na、nb为A、B阶次

L=400; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
xik=zeros(na,1); %白噪声初值ξ
etak=zeros(d+nb,1); %白噪声初值η
u=randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列ξ
eta=sqrt(0.25)*randn(L,1); %白噪声序列η

theta=[a(2:na+1);b]; %对象参数真值

thetae_1=zeros(na+nb+1,1); %参数估计初值
Rk_1=eye(na+nb+1);
for k=1:L
    phi=[-yk;uk(d:d+nb)];
    e(k)=a'*[xi(k);xik]-b'*etak(d:d+nb);
    y(k)=phi'*theta+e(k); %采集输出数据
    
    %随机牛顿算法
    R=Rk_1+(phi*phi'-Rk_1)/k;
    dR=det(R);
    if abs(dR)<10^(-6)  %避免矩阵R非奇异
        R=eye(na+nb+1);
    end
    IR=inv(R);
    thetae(:,k)=thetae_1+IR*phi*(y(k)-phi'*thetae_1)/k;
           
    %更新数据
    thetae_1=thetae(:,k);
    Rk_1=R;
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
        etak(i)=etak(i-1);
    end
    uk(1)=u(k);
    etak(1)=eta(k);  
    
    for i=na:-1:2
        yk(i)=yk(i-1);
        xik(i)=xik(i-1);
    end
    yk(1)=y(k);
    xik(1)=xi(k);
end
plot([1:L],thetae);
xlabel('k'); ylabel('参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);