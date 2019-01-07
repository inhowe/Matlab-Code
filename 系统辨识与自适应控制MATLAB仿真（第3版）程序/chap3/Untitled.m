%不需预选di的n阶SISO离散系统MRAS（用于参数估计）
%能被辨识的系统
clear all; 
close all;

% am=conv([1 0.8],[1 0.4])'
% bm=[3 2]'; %参考模型参数（参考模型中含有yr(k)，注意nb的使用！）
% am=[1 0.5 -0.6]';
% bm=[3 2]'; %参考模型参数（参考模型中含有yr(k)，注意nb的使用！）
am=conv([1 -0.2],[1 0.4])';
am=-am(2:3)
bm=[3]'; %参考模型参数（参考模型中含有yr(k)，注意nb的使用！）
thetam=[am;bm]; %参考模型参数向量
na=length(am); nb=length(bm); %可调参数个数

L=1000; %仿真长度
ypk=zeros(na,1); ymk=zeros(na,1); yrk=zeros(nb-1,1); epsilonk=zeros(na,1); %初值：ypk(i)表示yp(k-i)
yr=rands(L,1); %参考输入

G_1=eye(2*na+nb); lambda=1; %正定对称矩阵G(0)、λ

thetape_1=zeros(2*na+nb,1); %可调参数初值θpe(0)
for k=1:L
    time(k)=k;
    xm_1=[ymk;yr(k);yrk]; %构造参考模型数据向量
    ym(k)=thetam'*xm_1; %采集参考模型输出
    
    xpe_1=[ypk;yr(k);yrk;epsilonk]; %构造对象数据向量
    v0(k)=ym(k)-thetape_1'*xpe_1; %计算v0(k)
 
    G=G_1-G_1*xpe_1*xpe_1'*G_1/lambda/(1+xpe_1'*G_1*xpe_1/lambda); %计算G
    thetape(:,k)=thetape_1+G_1*xpe_1*v0(k)/(1+xpe_1'*G_1*xpe_1); %计算θpe
    
    yp(k)=thetape(1:na+nb,k)'*xpe_1(1:na+nb); %利用最新可调参数计算yp
    epsilon=ym(k)-yp(k); %广义输出误差ε
    
    %更新数据
    thetape_1=thetape(:,k);
    G_1=G;
    
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
plot(time,thetape(1:na,:));
xlabel('k'); ylabel('可调系统参数ap');
legend('a_p_1','a_p_2','a_p_3'); 
% axis([0 L -1 1.5]);
subplot(2,1,2);
plot(time,thetape(na+1:na+nb,:));
xlabel('k'); ylabel('可调系统参数bp');
legend('b_p_0','b_p_1'); 
% axis([0 L 0 4]);