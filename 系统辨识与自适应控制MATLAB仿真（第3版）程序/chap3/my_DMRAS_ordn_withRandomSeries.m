%随机序列下的不需预选di的n阶SISO离散系统MRAS（用于参数估计）
%被辨识系统格式：'_'表示下标
%              y(k)=Σam_i*y(k-i)+Σbm_i*u(k-i)
%输入信号必须平稳随机信号，阶跃信号不可用,因为阶跃信号无法充分激励系统的响应。
clf;
clear all;

EndTime=100000;%仿真截止时间
sampleTime=0.01;%采样间隔
t=0:sampleTime:EndTime;
L=EndTime/sampleTime; %仿真长度

num=[0.0005 0];
den=[1 -0.99];
% G=tf([0.98407],[13.407 1]);%不好辨识，离散化后分子系数太小，需要很长的辨识周期。
% Gz=c2d(G,sampleTime,'i')
G=tf(num,den,sampleTime);%测试系统
Gz=G
% G=tf([3 2 0 0],[1 -1 0.7 0.5],sampleTime);%例程系统，书本上的
% Gz=G
% seqGz=step(Gz,t);%Gz阶跃响应序列，不应该用
rng(1);%设置随机数种子，保证随机数可重复
yk=0;yk1=0;yk2=0;yk3=0;
% uk=idinput(L);%官方的生成-1/1随机序列
uk=fix(rands(L,1)+1);%生成01随机序列
% uk=sign(rands(L,1));%生成-1/1随机序列，01序列似乎更多情况下容易获得正确结果
uk1=0;uk2=0;uk3=0;
%生成被辨识系统的输出信号
for k=1:L
%     yk(k)=[1 -0.7 0.5]*[yk1 yk2 yk3]'+[3 2 0 0]*[uk(k) uk1 uk2 uk3]';
    %记得要改为和假定被测系统一致
    yk(k)=[den(2) 0 0]*[yk1 yk2 yk3]'+[num(1) num(2) 0]*[uk(k) uk1 uk2]';

    yk3=yk2;
    yk2=yk1;
    yk1=yk(k);
    uk3=uk2;
    uk2=uk1;
    uk1=uk(k);
end

na=length(den)-1; nb=length(num); %期望的可调参数个数：即na表示向量am的维数

ypk=zeros(na,1); ymk=zeros(na,1); yrk=zeros(nb-1,1); epsilonk=zeros(na,1); %初值：ypk(i)表示yp(k-i)
% yr=ones(L,1); %参考输入
yr=uk;%随机输入信号更容易辨识成功也更准确！

G_1=eye(2*na+nb); lambda=1; %正定对称矩阵G(0)、λ

thetape_1=zeros(2*na+nb,1); %可调参数初值θpe(0)
for k=1:L
    time(k)=k;
    ym(k)=yk(k); %采集参考模型输出，随机输入信号更容易辨识成功也更准确！
%     ym(k)=seqGz(k);
    
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
% legend('a_p_1','a_p_2','a_p_3'); 
% axis([0 L -1 1.5]);
subplot(2,1,2);
plot(time,thetape(na+1:na+nb,:));
xlabel('k'); ylabel('可调系统参数bp');
% legend('b_p_0','b_p_1'); 
% axis([0 L 0 4]);