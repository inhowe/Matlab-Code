function [E,F,G]=multidiophantine(a,b,c,N)
%********************************************************
  %功能：多步Diophanine方程的求解
  %调用格式：[E,F,G]=sindiophantine(a,b,c,N)（注：d=1）
  %输入参数：多项式A、B、C系数向量及预测步数（共4个）
  %输出参数：Diophanine方程的解E、F、G（共3个）
%********************************************************
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %A、B、C的阶次

%E、F、G的初值
E=zeros(N); E(1,1)=1; F(1,:)=conv(b,E(1,:)); 
if na>=nc
    G(1,:)=[c(2:nc+1) zeros(1,na-nc)]-a(2:na+1); %令c(nc+2)=c(nc+3)=...=0
else
    G(1,:)=c(2:nc+1)-[a(2:na+1) zeros(1,nc-na)]; %令a(na+2)=a(na+3)=...=0
end

%求E、G、F
for j=2:N
    for i=1:j-1
        E(j,i)=E(j-1,i);
    end
    E(j,j)=G(j-1,1);
    for i=2:na
        G(j,i-1)=G(j-1,i)-G(j-1,1)*a(i);
    end
    G(j,na)=-G(j-1,1)*a(na+1);
    F(j,:)=conv(b,E(j,:));
end