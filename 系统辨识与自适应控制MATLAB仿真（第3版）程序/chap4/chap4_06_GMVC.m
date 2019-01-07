%广义最小方差控制（显式控制律）
clear all; close all;

a=[1 -1.7 0.7]; b=[1 2]; c=[1 0.2]; d=4; %对象参数（无稳态误差）
%a=[1 -2 0.7]; b=[1 2]; c=[1 0.2]; d=4; %对象参数（有稳态误差）
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为多项式A、B、C阶次
nf=nb+d-1; ng=na-1; %nf、ng为多项式F、G的阶次

P=1; R=1; Q=2; %加权多项式
np=length(P)-1; nr=length(R)-1; nq=length(Q)-1;

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i);
yk=zeros(na,1); %输出初值
yrk=zeros(nc,1); %期望输出初值
xik=zeros(nc,1); %白噪声初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出
xi=sqrt(0.1)*randn(L,1); %白噪声序列

[e,f,g]=sindiophantine(a,b,c,d); %求解单步Diophantine方程
CQ=conv(c,Q); FP=conv(f,P); CR=conv(c,R); GP=conv(g,P); %CQ=C*Q
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb)+c*[xi(k);xik]; %采集输出数据
           
    u(k)=(-Q(1)*CQ(2:nc+nq+1)*uk(1:nc+nq)/b(1)-FP(2:np+nf+1)*uk(1:np+nf)...
        +CR*[yr(k+d:-1:k+d-min(d,nr+nc)); yrk(1:nr+nc-d)]...
        -GP*[y(k); yk(1:np+ng)])/(Q(1)*CQ(1)/b(1)+FP(1)); %求控制量
    
    %更新数据
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
    end
    if nc>0
        yrk(1)=yr(k);
        xik(1)=xi(k);
    end
end
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)');
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)');