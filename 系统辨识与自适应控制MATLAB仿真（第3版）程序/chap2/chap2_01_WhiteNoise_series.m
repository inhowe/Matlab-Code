%白噪声及有色噪声序列的产生
clear all; close all;

L=500; %仿真长度
d=[1 -1.5 0.7 0.1]; c=[1 0.5 0.2]; %D、C多项式的系数(可用roots命令求其根)
nd=length(d)-1; nc=length(c)-1; %nd、nc为D、C的阶次
xik=zeros(nc,1); %白噪声初值，相当于ξ(k-1)...ξ(k-nc)
ek=zeros(nd,1); %有色噪声初值
xi=randn(L,1); %randn产生均值为0，方差为1的高斯随机序列（白噪声序列）

for k=1:L
    e(k)=-d(2:nd+1)*ek+c*[xi(k);xik]; %产生有色噪声
    
    %数据更新
    for i=nd:-1:2
        ek(i)=ek(i-1);
    end
    ek(1)=e(k);
    
    for i=nc:-1:2
        xik(i)=xik(i-1);
    end
    xik(1)=xi(k);
end
subplot(2,1,1);
plot(xi);
xlabel('k'); ylabel('噪声幅值'); title('白噪声序列');
subplot(2,1,2);
plot(e);
xlabel('k'); ylabel('噪声幅值'); title('有色噪声序列');