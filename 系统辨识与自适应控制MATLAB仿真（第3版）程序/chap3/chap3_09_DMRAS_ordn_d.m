%需预选di的n阶SISO离散系统MRAS（用于参数估计）
clear all; close all;

am=[1 -0.7 0.5]'; bm=[3 2]'; %参考模型参数（参考模型中含有yr(k)，注意nb的使用！）
thetam=[am;bm]; %参考模型参数向量
na=length(am); nb=length(bm); %可调参数个数

L=400; %仿真长度
ypk=zeros(na,1); ymk=zeros(na,1); yrk=zeros(nb-1,1); epsilonk=zeros(na,1); %初值：ypk(i)表示yp(k-i)
yr=rands(L,1); %参考输入

G=1*eye(na+nb); G1=0.1*G; %注意G1的选择(此例选为常数矩阵)
D=-am+0.1; %系数di

thetapr_1=zeros(na+nb,1); %可调参数初值θpr
for k=1:L
    time(k)=k;
    xm_1=[ymk;yr(k);yrk]; %构造参考模型数据向量
    ym(k)=thetam'*xm_1; %采集参考模型输出
    
    xp_1=[ypk;yr(k);yrk]; %构造对象数据向量
    v0(k)=ym(k)-thetapr_1'*xp_1+D'*epsilonk; %计算v0(k)
    
    thetapp(:,k)=G1*xp_1*v0(k)/(1+xp_1'*(G+G1)*xp_1); %计算θpp
    thetapr(:,k)=thetapr_1+G*xp_1*v0(k)/(1+xp_1'*(G+G1)*xp_1); %计算θpr
    thetap(:,k)=thetapr(:,k)+thetapp(:,k); %θp
    
    yp(k)=thetap(:,k)'*xp_1; %利用最新可调参数计算yp
    epsilon=ym(k)-yp(k); %广义输出误差ε
    
    %更新数据
    thetapr_1=thetapr(:,k);
    
    for i=nb-1:-1:2  %注意参考模型中有yr(k)
        yrk(i)=yrk(i-1);
    end
    if nb>1
        yrk(1)=yr(k);
    end
    
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
plot(time,thetap(1:na,:));
xlabel('k'); ylabel('可调系统参数ap');
legend('a_p_1','a_p_2','a_p_3'); axis([0 L -1 1.5]);
subplot(2,1,2);
plot(time,thetap(na+1:na+nb,:));
xlabel('k'); ylabel('可调系统参数bp');
legend('b_p_0','b_p_1'); axis([0 L 0 4]);