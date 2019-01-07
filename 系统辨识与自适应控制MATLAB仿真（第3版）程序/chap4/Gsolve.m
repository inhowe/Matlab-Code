function G=Gsolve(a,b,d,N)
%****************************************************************
  %功能：改进GPC控制矩阵G的求解（CARIMA模型）
  %调用格式：G=Gsolve(a,b,d,N)
  %输入参数：多项式A、B（行向量）、纯滞后和预测长度（共4个）
  %输出参数：控制矩阵G
%****************************************************************
na=length(a)-1; nb=length(b)-1; %na是多项式A的阶次，nb是多项式B的阶次
a1=a(2:na+1);

G=zeros(N-d+1);
G(1,1)=b(1);
for j=2:N-d+1
    ab=0;
    for i=1:min(j-1,na)
        ab=ab+a1(i)*G(j-i,1);
    end
    if j<=nb+1
        b1j=b(j);
    else
        b1j=0;
    end
    G(j,1)=b1j-ab;
    for i=2:j
        G(j,i)=G(j-1,i-1);
    end
end