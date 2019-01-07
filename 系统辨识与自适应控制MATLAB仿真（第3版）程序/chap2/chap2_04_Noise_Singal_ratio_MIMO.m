%MIMO系统噪信比的计算
%“A(z)z(k)=B(z)u(k)+A(z)w(k)”（对应多变量辨识仿真实例2.12）
clear all;

a=[1 -0.8 -0.2 0.6]'; %多项式A
b(:,1,1)=[0 3 -3.5 -1.5 0]; %B11
b(:,1,2)=[0 1 -0.2 -0.5 0]; %B12
b(:,2,1)=[0 0 -4 -2 -1]; %B21
b(:,2,2)=[0 1 -1.5 0.5 0.2]; %B22（所有Bij维数需一致）
na=length(a)-1; Sz=size(b); nb=Sz(1)-1; %na、nb为A、B阶次
m=Sz(2); r=Sz(3); %m、r为系统输出输入维数

n=max(na,nb);
a0=[a;zeros(n-na,1)]; %高次项补0
for j=1:r
    b0(:,:,j)=[b(:,:,j); zeros(n-nb,m)];
end
deltau2=[1;1]; deltav2=0.1*[1;1]; %输入、白噪声方差

for i=1:n+1  %计算p、q的初值
    p(i,n+1)=a0(i);
    ii=0;
    for j=1:m
        for k=1:r
            ii=ii+1;
            q(i,n+1,ii)=b0(i,j,k);
        end
    end
end

for k=n:-1:1  %计算p、q
    for i=1:k
        p(i,k)=(p(1,k+1)*p(i,k+1)-p(k+1,k+1)*p(k+2-i,k+1))/p(1,k+1);
        for j=1:m*r
            q(i,k,j)=(p(1,k+1)*q(i,k+1,j)-q(k+1,k+1,j)*p(k+2-i,k+1))/p(1,k+1);
        end
    end
end

deltax2=zeros(m,1); ii=0;
for i=1:m  %求输入响应x的方差
    for j=1:r
        ii=ii+1;
        for k=1:n+1
            deltax2(i)=deltax2(i)+q(k,k,ii)^2*deltau2(j)/p(1,k);
        end
    end
end
deltax2=deltax2/a(1);
deltae2=deltav2;
ns=sqrt(deltae2./deltax2)  %噪信比（结果见MATLAB命令窗口）