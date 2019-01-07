%极点配置控制（PPC）（连续系统离散化）
clear all; close all;

%被控对象离散化
den=[1 1 0]; num=[1]; Ts=0.5; Td=0; %连续系统对象参数
sys=tf(num,den,'inputdelay',Td); %连续系统传递函数
dsys=c2d(sys,Ts,'zoh'); %离散化
[dnum,a]=tfdata(dsys,'v'); %提取离散系统数据
na=length(a)-1; b=dnum(2:na+1); nb=length(b)-1;
d=Td/Ts+1; %纯延时

%期望特性离散化
den=[1 2*0.7*1 1^2]; num=[1];
sys=tf(num,den);
dsys=c2d(sys,Ts,'zoh');
[dnum,Am]=tfdata(dsys,'v'); %提取Am
nam=length(Am)-1; %期望特征多项式阶次

%多项式B的分解
br=roots(b); %求B的根
b0=b(1); b1=1; %b0为B-；b1为B+
Val=0.8; %通过修改临界值，确定B零点是否对消（零点绝对值小于临界值，则被抵消）
for i=1:nb %分解B-、B+
    if abs(br(i))>=Val
        b0=conv(b0,[1 -br(i)]);
    else
        b1=conv(b1,[1 -br(i)]);
    end
end

Bm1=sum(Am)/sum(b0);Bm=Bm1*b0; %确定多项式Bm

%确定多项式A0
na0=2*na-1-nam-(length(b1)-1); %观测器最低阶次
A0=1;
for i=1:na0
    A0=conv(A0,[1 0.5]); %生成观测器
end

%计算Diophantine方程，得到F、G、R
[F1,G]=diophantine(a,b0,d,A0,Am); %注意，此处为b0
F=conv(F1,b1); R=Bm1*A0; 
nf=length(F)-1; ng=length(G)-1; nr=length(R)-1;

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
yrk=zeros(na,1); %期望输出初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出

for k=1:L
    time(k)=k*Ts;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb); %采集输出数据
    
    u(k)=(-F(2:nf+1)*uk(1:nf)+R*[yr(k+d:-1:k+d-min(d,nr));yrk(1:nr-d)]-G*[y(k);yk(1:ng)])/F(1); %求控制量

    %更新数据
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
        yrk(i)=yrk(i-1);
    end
    yk(1)=y(k);
    yrk(1)=yr(k);
end
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('t'); ylabel('y_r(t)、y(t)');
legend('y_r(t)','y(t)');
subplot(2,1,2);
plot(time,u);
xlabel('t'); ylabel('u(t)');