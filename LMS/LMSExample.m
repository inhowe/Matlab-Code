%function main()
close  all

% 周期信号的产生 
t=0:0.1:99;
xs=10*sin(0.5*t);
figure;
subplot(2,1,1);
plot(t,xs);grid;
ylabel('幅值');
title('it{输入周期性信号}');

% 噪声信号的产生
randn('state',sum(100*clock));
xn=randn(1,size(t,2));
subplot(2,1,2);
plot(t,xn);grid;
ylabel('幅值');
xlabel('时间');
title('it{随机噪声信号}');

% 信号滤波
xn = xs+xn;
xn = xn.' ;   % 输入信号序列
dn = xs.' ;   % 预期结果序列
M  = 10   ;   % 滤波器的阶数

rho_max = max(eig(xn*xn.'));   % 输入信号相关矩阵的最大特征值
mu = rand()*(1/rho_max)   ;    % 收敛因子 0 < mu < 1/rho

[yn,W,en] = LMS(xn,dn,M,mu);

% 绘制滤波器输入信号
figure;
subplot(2,1,1);
plot(t,xn);grid;
ylabel('幅值');
xlabel('时间');
title('it{滤波器输入信号}');

% 绘制自适应滤波器输出信号
subplot(2,1,2);
plot(t,yn);grid;
ylabel('幅值');
xlabel('时间');
title('it{自适应滤波器输出信号}');

% 绘制自适应滤波器输出信号,预期输出信号和两者的误差
figure 
plot(t,yn,'b',t,dn,'g',t,dn-yn,'r');grid;
axis([0 100 -12 12]);
legend('自适应滤波器输出','预期输出','误差');
ylabel('幅值');
xlabel('时间');
Str=num2str(M);
title(['it{自适应滤波器}' ' 阶数M=' Str]);

%--------------------------------------------------
% mu = rand()*(1/rho_max);    % 换一个收敛因子重新计算以进行比较
M  = 20;                    % 换一个滤波器的阶数重新计算以进行比较
[yn,W,en] = LMS(xn,dn,M,mu);

figure 
plot(t,yn,'b',t,dn,'g',t,dn-yn,'r');grid;
axis([0 100 -12 12]);
legend('自适应滤波器输出','预期输出','误差');
ylabel('幅值');
xlabel('时间');
Str=num2str(M);
title(['it{自适应滤波器}' ' 阶数M=' Str]);
% title(['it{自适应滤波器}' ' mu=' Str_mu]);