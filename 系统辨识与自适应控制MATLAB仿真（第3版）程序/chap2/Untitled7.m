clear all; close all;
sampleTime=0.01;%采样时间
simTime=10;%仿真时间
% 假定一个传函为待辨识系统
num=[0.9];
den=[1 2];
G=tf(num,den)
Gz=c2d(G,sampleTime,'i') % 离散化得到脉冲传函，需要说明的是，有时候分子会出现很小的值是正常现象，可以忽略为0，因此该项的辨识也是有可能不精确的
% Gz=tf(num,den,sampleTime)% 直接假定一个离散化的传函，稳定的系统才可正确辨识（单位圆内）
%此时可以得到差分表达式：y/u=Gz
%化为：A(z)y(k)=B(z)u(k)*z^-d
%其中：A(z)=1+a1*z^-1+a2*z^-2+...
yt=step(Gz,0:sampleTime:simTime);%得到阶跃输出序列
step(Gz)
figure
step(feedback(Gz,1))
[num_Gz,den_Gz]=tfdata(Gz,'v');

%确定性系统的递推梯度校正参数估计（RGC）
d=1; %对象参数
n_den_Gz=length(den_Gz)-1; n_num_Gz=length(num_Gz)-1-d+1; %na、nb为A、B阶次

L=simTime/sampleTime; %仿真长度
uk=zeros(d+n_num_Gz,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(n_den_Gz,1); %输出初值
u=ones(L,1); %输入采用白噪声序列
% xi=sqrt(0.1)*randn(L,1);%噪声
% theta=[den_Gz(2:na+1);num_Gz]; %对象参数真值

thetae_1=zeros(n_den_Gz+n_num_Gz+1,1); %参数估计初值
alpha=1; %范围(0,2)
c=0.1; %修正因子
for k=1:L
    phi=[-yk;uk(d:d+n_num_Gz)];
    y(k)=yt(k); %采集输出数据
    
    thetae(:,k)=thetae_1+alpha*phi*(y(k)-phi'*thetae_1)/(phi'*phi+c); %递推梯度校正算法
    
    %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+n_num_Gz:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=n_den_Gz:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
thetae_1'
% plot([1:L],thetae);
% xlabel('k'); ylabel('参数估计a、b');
% legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);