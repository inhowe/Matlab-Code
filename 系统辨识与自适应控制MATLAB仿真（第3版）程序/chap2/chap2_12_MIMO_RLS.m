%MIMO系统递推最小二乘参数估计（本程序针对2入2出系统）
clear all; close all;

a=[1 -0.8 -0.2 0.6]'; %多项式A
b(:,1,1)=[3 -3.5 -1.5 0]; %Bij（维数需相同）
b(:,1,2)=[1 -0.2 -0.5 0];
b(:,2,1)=[0 -4 -2 -1];
b(:,2,2)=[1 -1.5 0.5 0.2];
d=[1 1;2 1]; %Bij每一项的滞后d
nb=[2 2;2 3]; %Bij每一项的阶次
na=length(a)-1; Sz=size(b); r=Sz(3); m=Sz(2); %r、m分别为系统输入输出维数
mb=max(max(nb+d)); %用于数据更新

L=400; %仿真长度
uk=zeros(mb,r); %输入初值
yk=zeros(na,m); %输出初值
xik=zeros(na,m); %白噪声初值
xi(:,1)=sqrt(0.1)*randn(L,1); xi(:,2)=sqrt(0.1)*randn(L,1); %白噪声
u(:,1)=randn(L,1); u(:,2)=randn(L,1); %输入量

%构造向量θ（对象参数真值）
theta1=[]; theta2=[];
for j=1:r
    theta1=[theta1;b(d(1,j):d(1,j)+nb(1,j),1,j)];
    theta2=[theta2;b(d(2,j):d(2,j)+nb(2,j),2,j)];
end
theta=[a(2:na+1)' theta1' theta2']';

thetae_1=zeros(na+sum(sum(nb))+4,1); %参数估计初值
P=10^6*eye(na+sum(sum(nb))+4);
for k=1:L
    %组建Φ(k)
    u1k=[]; u2k=[];
    for j=1:r
        u1k=[u1k;uk(d(1,j):d(1,j)+nb(1,j),j)];
        u2k=[u2k;uk(d(2,j):d(2,j)+nb(2,j),j)];
    end
    phi=[-yk', [u1k';zeros(size(u1k'))], [zeros(size(u2k'));u2k']];
    
    e=[xi(k,:)',xik']*a;
    y(:,k)=phi*theta+e; %采集输出数据
    
    %递推最小二乘法
    K=P*phi'*inv(phi*P*phi'+eye(m));
    thetae(:,k)=thetae_1+K*(y(:,k)-phi*thetae_1);
    P=(eye(na+sum(sum(nb))+4)-K*phi)*P;
        
    %更新数据
    thetae_1=thetae(:,k);
    
    for j=1:r
        for i=mb:-1:2
            uk(i,j)=uk(i-1,j);
        end
        uk(1,j)=u(k,j);
    end
    
    for j=1:m
        for i=na:-1:2
            yk(i,j)=yk(i-1,j);
            xik(i,j)=xik(i-1,j);
        end
        yk(1,j)=y(j,k);
        xik(1,j)=xi(k,j);
    end
end
figure(1)
plot([1:L],thetae(1:na,:));
xlabel('k'); ylabel('参数估计a');
legend('a(1)','a(2)','a(3)'); axis([0 L -1 1]);
figure(2)
plot([1:L],thetae(na+1:na+nb(1,1)+1,:));
xlabel('k'); ylabel('参数估计b11');
legend('b_1_1(0)','b_1_1(1)','b_1_1(2)'); axis([0 L -4 4]);
figure(3)
plot([1:L],thetae(na+nb(1,1)+2:na+sum(nb(1,:))+2,:));
xlabel('k'); ylabel('参数估计b12');
legend('b_1_2(0)','b_1_2(1)','b_1_2(2)'); axis([0 L -1 2]);
figure(4)
plot([1:L],thetae(na+sum(nb(1,:))+3:na+sum(nb(1,:))+nb(2,1)+3,:));
xlabel('k'); ylabel('参数估计b21');
legend('b_2_1(0)','b_2_1(1)','b_2_1(2)'); axis([0 L -4.5 0]);
figure(5)
plot([1:L],thetae(na+sum(nb(1,:))+nb(2,1)+4:na+sum(sum(nb))+4,:));
xlabel('k'); ylabel('参数估计b22');
legend('b_2_2(0)','b_2_2(1)','b_2_2(2)','b_2_2(3)'); axis([0 L -2 2]);