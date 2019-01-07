%基于CARMA模型的JGPC（C!=1）（对象参数未知）
clear all; close all;

a=[1 -2 1.1]; b=[1 2]; c=[1 -2.5]; d=4; %对象参数（C可以不稳定）
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为多项式A、B、C阶次

N=8; gamma=1*eye(N-d+1); alpha=0.7; %预测长度、控制加权矩阵、输出柔化系数

L=600; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
xik=zeros(nc,1); %白噪声初值
xiek=zeros(nc,1); %白噪声估计初值
w=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %设置值
xi=sqrt(0.01)*randn(L,1); %白噪声序列

%RELS初值
thetae_1=0.001*ones(na+nb+1+nc,1);
P=10^6*eye(na+nb+1+nc);
lambda=1; %遗忘因子[0.9 1]
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb)+c*[xi(k);xik]; %采集输出数据
    
    %递推增广最小二乘法
    phie=[-yk(1:na);uk(d:d+nb);xiek];
    K=P*phie/(lambda+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1+nc)-K*phie')*P/lambda;
    
    xie=y(k)-phie'*thetae(:,k);%白噪声的估计值
    
    %提取辨识参数
    ae=[1 thetae(1:na,k)']; be=thetae(na+1:na+nb+1,k)'; ce=[1 thetae(na+nb+2:na+nb+1+nc,k)'];
    
    G=Gsolve(ae,be,d,N); %求解控制矩阵G
    
    %计算模型预测输出向量Ym
    ym(k)=y(k);
    for j=1:N
        ym(k+j)=-ae(2:na+1)*[ym(k+j-1:-1:k+j-min(j,na))';yk(1:na-j)];
        for i=0:nb
            if j-d-i>=0
                ym(k+j)=ym(k+j)+be(i+1)*uk(1);
            else
                ym(k+j)=ym(k+j)+be(i+1)*uk(d+i-j);
            end
        end
        for i=0:nc
            if j-i==0
                ym(k+j)=ym(k+j)+ce(i+1)*xi(k);
            elseif j-i<0
                ym(k+j)=ym(k+j)+ce(i+1)*xik(i-j);
            end
        end
    end
    Ym=[ym(k+d:k+N)]'; %构造向量Ym
    
    %参考轨迹向量Yr
    yr(k+d-1)=ym(k+d-1);
    for i=0:N-d
        yr(k+d+i)=alpha*yr(k+d+i-1)+(1-alpha)*w(k+d);
    end
    Yr=[yr(k+d:k+N)]'; %构造向量Yr
    
    %求控制律
    dU=inv(G'*G+gamma)*G'*(Yr-Ym); %ΔU
    du(k)=dU(1); u(k)=uk(1)+du(k);
    
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
    if nc>0
        xik(1)=xi(k);
        xiek(1)=xie;
    end
end
figure(1);
subplot(2,1,1);
plot(time,w(1:L),'r:',time,y);
xlabel('k'); ylabel('w(k)、y(k)');
legend('w(k)','y(k)'); axis([0 L -20 20]);
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -4 4]);
figure(2)
plot([1:L],thetae);
xlabel('k'); ylabel('辨识参数a、b、c');
legend('a_1','a_2','b_0','b_1','c_1'); axis([0 L -3 3]);