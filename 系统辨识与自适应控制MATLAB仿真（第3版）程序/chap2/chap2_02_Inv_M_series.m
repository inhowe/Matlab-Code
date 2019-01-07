%M序列及逆M序列的产生
clear all; close all;

L=60; %M序列长度
x1=1; x2=1; x3=1; x4=0; %移位寄存器初值xi-1、xi-2、xi-3、xi-4
S=1; %方波初值

for k=1:L
    M(k)=xor(x3,x4); %进行异或运算,产生M序列
    IM=xor(M(k),S); %进行异或运算,产生逆M序列
    if IM==0
        u(k)=-1;
    else
        u(k)=1;
    end
    
    S=not(S); %产生方波
    
    x4=x3; x3=x2; x2=x1; x1=M(k); %寄存器移位
    
end
subplot(2,1,1);
stairs(M); grid;
axis([0 L/2 -0.5 1.5]); xlabel('k'); ylabel('M序列幅值'); title('M序列');
subplot(2,1,2);
stairs(u); grid;
axis([0 L -1.5 1.5]); xlabel('k'); ylabel('逆M序列幅值'); title('逆M序列');