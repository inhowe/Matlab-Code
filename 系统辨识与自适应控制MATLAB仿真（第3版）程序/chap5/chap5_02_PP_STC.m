%极点配置间接自校正控制
clear all; close all;

a=[1 -1.6065 0.6065]; b=[0.1065 0.0902]; c=[1 0.5]; d=10; %对象参数
Am=[1 -1.3205 0.4966]; %期望闭环特征多项式
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为多项式A、B、C
nam=length(Am)-1; %Am阶次
nf=nb+d-1; ng=na-1;

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
yrk=zeros(na,1); %期望输出初值
xik=zeros(nc,1); %白噪声初值
xiek=zeros(nc,1); %白噪声估计初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出
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
    
    xie=y(k)-phie'*thetae(:,k); %白噪声的估计值
    
    %提取辨识参数
    ae=[1 thetae(1:na,k)']; be=thetae(na+1:na+nb+1,k)'; ce=[1 thetae(na+nb+2:na+nb+1+nc,k)'];
    if nc>0
        if abs(ce(2))>0.8
            ce(2)=sign(ce(2))*0.8;
        end
    end
    
    %多项式B的分解
    br=roots(be); %求B的根
    b0=be(1); b1=1; %b0为B-；b1为B+
    Val=0.9; %通过修改临界值，确定B零点是否对消（零点绝对值小于临界值，则被抵消）
    for i=1:nb %分解B-、B+
        if abs(br(i))>=Val
            b0=conv(b0,[1 -br(i)]);
        else
            b1=conv(b1,[1 -br(i)]);
        end
    end
    
    Bm1=sum(Am)/sum(b0); %确定多项式Bm'
    
    %确定多项式A0
    %A0=ce; %可取A0=C
    na0=2*na-1-nam-(length(b1)-1); %观测器最低阶次
    A0=1;
    for i=1:na0
        A0=conv(A0,[1 0.5]); %生成观测器
    end

    %计算Diophantine方程，得到F、G、R
    [F1,G]=diophantine(ae,b0,d,A0,Am);
    F=conv(F1,b1); R=Bm1*A0; 
    nr=length(R)-1;

    u(k)=(-F(2:nf+1)*uk(1:nf)+R*[yr(k+d:-1:k+d-min(d,nr));yrk(1:nr-d)]-G*[y(k);yk(1:ng)])/F(1);%求控制量
    
    %更新数据
    thetae_1=thetae(:,k);
    
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
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L -20 20]);
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -40 40]);
figure(2)
subplot(211)
plot([1:L],thetae(1:na,:),[1:L],thetae(na+nb+2:na+nb+1+nc,:));
xlabel('k'); ylabel('参数估计a、c');
legend('a_1','a_2','c_1'); axis([0 L -2 2]);
subplot(212)
plot([1:L],thetae(na+1:na+nb+1,:));
xlabel('k'); ylabel('参数估计b');
legend('b_0','b_1'); axis([0 L 0 0.15]);