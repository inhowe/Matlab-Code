%递推随机牛顿参数估计（RSNA）
clear all; close all;
sampleTime=0.01;%采样时间
simTime=10;%仿真时间
% 假定一个传函为待辨识系统
num=[1 0.5];
den=[1 -0.1 -0.12];
% G=tf(num,den)
% Gz=c2d(G,sampleTime,'i') % 离散化得到脉冲传函，需要说明的是，有时候分子会出现很小的值是正常现象，可以忽略为0，因此该项的辨识也是有可能不精确的
Gz=tf(num,den,sampleTime)% 直接假定一个离散化的传函，稳定的系统才可正确辨识（单位圆内）
%此时可以得到差分表达式：y/u=Gz
%化为：A(z)y(k)=B(z)u(k)*z^-d
%其中：A(z)=1+a1*z^-1+a2*z^-2+...
yt=step(Gz,0:sampleTime:simTime);%得到阶跃输出序列
step(Gz)
figure
step(feedback(Gz,1))
figure
[num_Gz,den_Gz]=tfdata(Gz,'v');

% a=[1 -0.1 -0.12]'; b=[1 0.5]'; 
d=1; %对象参数
n_Gz_num=length(num_Gz)-1-d;
n_Gz_den=length(den_Gz)-1;  %na、nb为A、B阶次
L=simTime/sampleTime; %仿真长度
u=ones(L,1); %输入采用白噪声序列
uk=zeros(d+n_Gz_num,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(n_Gz_den,1); %输出初值
thetae_1=zeros(n_Gz_den+n_Gz_num+1,1); %参数估计初值
Rk_1=eye(n_Gz_den+n_Gz_num+1);
for k=1:L
    phi=[-yk;uk(d:d+n_Gz_num)];
    y(k)=yt(k); %采集输出数据
    
    %随机牛顿算法
    R=Rk_1+(phi*phi'-Rk_1)/k;
    dR=det(R);
    if abs(dR)<10^(-6)  %避免矩阵R非奇异
        R=eye(n_Gz_den+n_Gz_num+1);
    end
    IR=inv(R);
    thetae(:,k)=thetae_1+IR*phi*(y(k)-phi'*thetae_1)/k;
           
    %更新数据
    thetae_1=thetae(:,k);
    Rk_1=R;
    
    for i=d+n_Gz_num:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=n_Gz_den:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
plot([1:L],thetae);
xlabel('k'); ylabel('参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);