%自校正PID控制（二阶系统，对象参数未知）
clear all; close all;

a=[1 -1.6065 0.6065]; b=[0.1065 0.0902]; d=3; Am=[1 -1.3205 0.4966]; %对象参数及期望闭环特征多项式
na=length(a)-1; nb=length(b)-1; nam=length(Am)-1; %na、nb、nc、nam为多项式A、B、C、Am阶次
nf1=nb+d+2-(na+1)+1; ng=2; %nf1=nf+1

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4,1)]; %期望输出
e=2*ones(L,1); %常值干扰

%RLS初值
thetae_1=0.001*ones(na+nb+1,1);
P=10^6*eye(na+nb+1);
lambda=1; %遗忘因子[0.9 1]
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb)+e(k); %采集输出数据
    
    %递推最小二乘法
    phie=[-yk(1:na);uk(d:d+nb)];
    K=P*phie/(lambda+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1)-K*phie')*P/lambda;
    
    %提取辨识参数
    ae=[1 thetae(1:na,k)']; be=thetae(na+1:na+nb+1,k)';
    
    %计算Diophantine方程，得到F、G、R
    [F,G]=diophantine(conv(ae,[1 -1]),be,d,1,Am); %A0=1
    F1=conv(F,[1 -1]); R=sum(G);
    
    u(k)=(-F1(2:nf1+1)*uk(1:nf1)+R*yr(k)-G*[y(k);yk(1:ng)])/F1(1);%求控制量
    
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
figure(1);
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L -20 20]);
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -40 20]);
figure(2)
subplot(211)
plot([1:L],thetae(1:na,:));
xlabel('k'); ylabel('参数估计a');
legend('a_1','a_2'); axis([0 L -2 2]);
subplot(212)
plot([1:L],thetae(na+1:na+nb+1,:));
xlabel('k'); ylabel('参数估计b');
legend('b_0','b_1'); axis([0 L 0 0.15]);