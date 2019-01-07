%最小方差间接自校正控制
clear all; close all;

a=[1 -1.7 0.7]; b=[1 0.5]; c=[1 0.2]; d=4; %对象参数
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为多项式A、B、C阶次
nf=nb+d-1; %nf为多项式F的阶次

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i);
yk=zeros(na,1); %输出初值
yrk=zeros(nc,1); %期望输出初值
xik=zeros(nc,1); %白噪声初值
xiek=zeros(nc,1); %白噪声估计初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出
xi=sqrt(0.1)*randn(L,1); %白噪声序列

%RELS初值设置
thetae_1=0.001*ones(na+nb+1+nc,1); %非常小的正数（这里不能为0）
P=10^6*eye(na+nb+1+nc);
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb)+c*[xi(k);xik]; %采集输出数据
        
    %递推增广最小二乘法
    phie=[-yk;uk(d:d+nb);xiek];
    K=P*phie/(1+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1+nc)-K*phie')*P;
    
    xie=y(k)-phie'*thetae(:,k); %白噪声的估计值
    
    %提取辨识参数
    ae=[1 thetae(1:na,k)']; be=thetae(na+1:na+nb+1,k)'; ce=[1 thetae(na+nb+2:na+nb+1+nc,k)'];
    if abs(be(2))>0.9
        be(2)=sign(ce(2))*0.9; %MVC算法要求B稳定
    end
    if abs(ce(2))>0.9
        ce(2)=sign(ce(2))*0.9;
    end
    
    [e,f,g]=sindiophantine(ae,be,ce,d); %求解单步Diophantine方程        
    u(k)=(-f(2:nf+1)*uk(1:nf)+ce*[yr(k+d:-1:k+d-min(d,nc));yrk(1:nc-d)]-g*[y(k);yk(1:na-1)])/f(1); %求控制量
    
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
        yrk(i)=yrk(i-1);
        xik(i)=xik(i-1);
        xiek(i)=xiek(i-1);
    end
    if nc>0
        yrk(1)=yr(k);
        xik(1)=xi(k);
        xiek(1)=xie;
    end
end
figure(1);
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L -20 20]);
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -40 40]);
figure(2)
subplot(211)
plot([1:L],thetae(1:na,:));
xlabel('k'); ylabel('参数估计a');
legend('a_1','a_2'); axis([0 L -3 2]);
subplot(212)
plot([1:L],thetae(na+1:na+nb+1+nc,:));
xlabel('k'); ylabel('参数估计b、c');
legend('b_0','b_1','c_1'); axis([0 L 0 1.5]);