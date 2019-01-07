%二阶SISO离散系统MRAS（用于参数估计）
clear all; close all;

am=[1 -0.2]'; bm=[2]'; d=1; %参考模型参数
na=length(am); nb=length(bm)-1; %na、nb为A、B阶次

L=400; %仿真长度
ypk=zeros(na,1); %初值：ypk(i)表示yp(k-i)
ymk=zeros(na,1);
yrk=zeros(nb+d,1);
epsilonk=zeros(na,1);
yr=ones(L,1); %参考输入

ap_1=zeros(na,1); bp_1=zeros(nb+1,1); %可调参数初值

alpha=ones(na,1); beta=ones(nb+1,1); %α>0、β>0
D=zeros(na,1); D(1)=4*(1+am(2))/(3+am(1)+am(2)); D(2)=D(1)-1; %系数di
for k=1:L
    time(k)=k;
    ym(k)=am'*ymk+bm'*yrk(d:d+nb); %采集参考模型输出
    
    yp0(k)=ap_1'*ypk+bp_1'*yrk(d:d+nb);
    v0(k)=ym(k)-yp0(k)+D'*epsilonk; %计算v0(k)
    for i=1:na
        ap(i,k)=ap_1(i)+alpha(i)*ypk(i)*v0(k)/(1+alpha'*ypk.^2+beta*yrk(d:d+nb).^2); %计算ap(k)
    end
    for i=1:nb+1
        bp(i,k)=bp_1(i)+beta(i)*yrk(i)*v0(k)/(1+alpha'*ypk.^2+beta*yrk(d:d+nb).^2); %计算bp(k)
    end
    
    yp(k)=ap(:,k)'*ypk+bp(:,k)'*yrk(d:d+nb); %计算yp
    epsilon=ym(k)-yp(k); %广义输出误差ε
    
    %更新数据
    ap_1=ap(:,k); bp_1=bp(:,k);
    
    for i=d+nb:-1:2
        yrk(i)=yrk(i-1);
    end
    yrk(1)=yr(k);
    
    for i=na:-1:2
        ypk(i)=ypk(i-1);
        ymk(i)=ymk(i-1);
        epsilonk(i)=epsilonk(i-1);
    end
    ypk(1)=yp(k);
    ymk(1)=ym(k);
    epsilonk(1)=epsilon;
end
subplot(2,1,1);
plot(time,ym,'r',time,yp,':');
xlabel('k'); ylabel('y_m(k)、y_p(k)');
legend('y_m(k)','y_p(k)');
subplot(2,1,2);
plot(time,[ap;bp]);
xlabel('k'); ylabel('对象参数ap、bp');
legend('a_p_1','a_p_2','b_p_1');