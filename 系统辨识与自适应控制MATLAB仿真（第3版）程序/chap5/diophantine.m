function [F1,G]=diophantine(A,B,d,A0,Am)
%***********************************************************************
  %功能：Diophanine方程的求解
  %调用格式：[F1,G]=diophantine(A,B,d,A0,Am)
  %输入参数：多项式A、B系数向量、纯延迟d、多项式A0、Am系数向量（行向量）
  %输出参数：Diophanine方程的解F1、G（行向量）
%***********************************************************************  
dB=[zeros(1,d) B];
na=length(A)-1; nd=length(dB)-1;
T1=conv(A0,Am); nt=length(T1); T=[T1';zeros(na+nd-nt,1)];

%得到Sylvester 矩阵
AB=zeros(na+nd);
for i=1:na+1
    for j=1:nd
        AB(i+j-1,j)=A(i);
    end
end
for i=1:nd+1
    for j=1:na
        AB(i+j-1,j+nd)=dB(i);
    end
end
%得到F1,G
L=(AB)\T;
F1=[ L(1:nd)]';
G=[ L(nd+1:na+nd)]';