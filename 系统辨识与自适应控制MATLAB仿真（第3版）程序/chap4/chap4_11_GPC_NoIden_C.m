%Clarke广义预测控制（C!=1）（对象参数已知）
%N1=d、N、Nu取不同的值
clear all; close all;

a=[1 -2 1.1]; b=[1 2]; c=[1 0.5]; d=4; %对象参数
na=length(a)-1; b=[zeros(1,d-1) b]; nb=length(b)-1; nc=length(c)-1;%多项式A、B、C阶次（因d!=1，对b添0）
aa=conv(a,[1 -1]); naa=na+1;%aa的阶次

N1=d; N=8; Nu=2; %最小输出长度、预测长度、控制长度
gamma=1*eye(Nu); alpha=0.7; %控制加权矩阵、输出柔化系数

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
duk=zeros(d+nb,1); %控制增量初值
dufk=zeros(max(nb,nc),1); %滤波控制增量初值
yk=zeros(naa,1); %输出初值
yfk=zeros(max(na,nc),1); %滤波输出初值
xik=zeros(nc,1); %白噪声初值
w=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %设定值
xi=sqrt(0.01)*randn(L,1); %白噪声序列

%求解多步Diophantine方程并构建F1、F2、G
[E,F,G]=multidiophantine(aa,b,c,N);
G=G(N1:N,:);
F1=zeros(N-N1+1,Nu); F2=zeros(N-N1+1,nb);
for i=1:N-N1+1
    for j=1:min(i,Nu);     F1(i,j)=F(i+N1-1,i+N1-1-j+1);  end
    for j=1:nb;            F2(i,j)=F(i+N1-1,i+N1-1+j);    end
end

for k=1:L
    time(k)=k;
    y(k)=-aa(2:naa+1)*yk+b*duk(1:nb+1)+c*[xi(k);xik]; %采集输出数据
    yf(k)=-c(2:nc+1)*yfk(1:nc)+y(k);
    Yfk=[yf(k);yfk(1:na)]; %构建向量Yf(k)
    dUfk=dufk(1:nb); %构建向量ΔUf(k-j)
    
    %参考轨迹
    yr(k)=y(k);
    for i=1:N
        yr(k+i)=alpha*yr(k+i-1)+(1-alpha)*w(k+d);
    end
    Yr=[yr(k+N1:k+N)]'; %构建向量Yr(k)
           
    %求控制律
    dUf=inv(F1'*F1+gamma)*F1'*(Yr-F2*dUfk-G*Yfk); %ΔUf
    duf(k)=dUf(1);
    du(k)=c*[duf(k); dufk(1:nc)];
    u(k)=uk(1)+du(k);
    
    %更新数据
    for i=1+nb:-1:2
        uk(i)=uk(i-1);
        duk(i)=duk(i-1);
    end
    uk(1)=u(k);
    duk(1)=du(k);
    
    for i=max(nb,nc):-1:2
        dufk(i)=dufk(i-1);
    end
    dufk(1)=duf(k);
    
    for i=naa:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
    
    for i=max(na,nc):-1:2
        yfk(i)=yfk(i-1);
    end
    yfk(1)=yf(k);
    
    for i=nc:-1:2
        xik(i)=xik(i-1);
    end
    if nc>0
        xik(1)=xi(k);
    end
end
subplot(2,1,1);
plot(time,w(1:L),'r:',time,y);
xlabel('k'); ylabel('w(k)、y(k)');
legend('w(k)','y(k)');
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)');