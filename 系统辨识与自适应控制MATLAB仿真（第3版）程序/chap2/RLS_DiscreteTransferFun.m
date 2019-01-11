% 系统辨识传递函数
% 将传递函数化为离散形式（AR模型），然后辨识这个AR模型的参数
% 离散的系统，极点必须全部位于单位圆内，分母为首一式，才可以辨识
% 可以开环不稳定但必须闭环稳定？
clear all; close all;
sampleTime=0.01;%采样时间
simTime=10;%仿真时间
% 假定一个传函为待辨识系统
num=[0.05 -0.03];
den=[1 -0.5 -0.5];
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

%准备工作（前2个变量建议人工核对一下！）
d=1;%d>0,d为待辨识系统（z域）的分子中z的最小指数（一般为0）+1
n_Gz_num=length(num_Gz)-1-d;%n_Gz_num+d等于待辨识系统（z域）的分子阶数
n_Gz_den=length(den_Gz)-1;%n_Gz_den为待辨识系统（z域）的分母阶数
L=simTime/sampleTime; %仿真（辨识）长度
u=ones(L,1); %生成输入序列，因为一般用的是阶跃信号

%准备递推最小二乘参数估计（RLS）
uk=zeros(d+n_Gz_num,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(n_Gz_den,1); %输出初值
thetae_1=zeros(n_Gz_den+n_Gz_num+1,1); %thetae初值
P=10^6*eye(n_Gz_den+n_Gz_num+1); 
for k=1:L
    % 从“测量”数据中取出一个测量点
    phi=[-yk;uk(d:d+n_Gz_num)]; %此处phi为列向量
    y(k)=yt(k); %采集输出数据
   
    %递推最小二乘法
    K=P*phi/(1+phi'*P*phi);
    thetae(:,k)=thetae_1+K*(y(k)-phi'*thetae_1);%使用thetae(:,k)添加新的向量进数组，为了保留历史辨识结果便于绘图
    P=(eye(n_Gz_den+n_Gz_num+1)-K*phi')*P;
    
    %更新数据
    thetae_1=thetae(:,k);
    
    %剔除旧的“测量”点，类似滑动窗口
    for i=d+n_Gz_num:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=n_Gz_den:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
disp('最新的估计值为'),thetae_1
hold on
tmp=size(thetae);%获取有多少个参数
for i=1:tmp(1)
    plot(thetae(i,:))
end
title('估计参数：分别为a0-an，b0-bn')