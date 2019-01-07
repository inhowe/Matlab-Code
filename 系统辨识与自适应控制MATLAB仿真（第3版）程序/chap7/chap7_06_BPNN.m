% BP神经网络辨识
clear all; close all;

[t,x] = ode45(@nonsys,[0 20],[1.0 1.0]); % 解微分方程
y = sin(x(:,1)+x(:,2));
[N,n] = size(x); % N为数据个数，n为系统输入的维数

% 设置BP网络参数
m = 10; % m为隐含层节点数
eta = 0.5; % 学习速率
alpha = 0.05; % 动量因子
E = 0.02; % 全局误差精度
% w1k1 = rands(m,n); % 输入层至隐含层权值的初值: w1ki表示w1(k-i)
w1k1 = [-0.7027    0.9150
    0.3162   -0.9486
    0.2680    0.9422
   -0.5414   -0.4048
   -0.6355    0.0501
   -0.6673    0.7247
   -0.7008    0.7928
   -0.5945   -0.6220
    0.9099    0.3214
   -0.9682    0.8825];
w1k2 = w1k1; 
% w2k1 = rands(1,m); % 隐含层至输出层权值的初值: w2ki表示w2(k-i)
w2k1 = [0.9514 -0.7841 -0.6422 0.4931 -0.9011 -0.8574 -0.0217 0.6998 0.9941 -0.9912];
w2k2 = w2k1; 

eg1 = 100*E; 
eg = 10; % 初始化全局误差
num = 0; % 初始化训练步数
M = 100; % 最大训练次数
tic
while(eg >= E) % 由全局误差控制训练次数
%while(num < M) % 直接设定训练次数
    num = num+1;
    es(num) = 0;
    for k = 1:N        
        % 计算BP网络输出
        O1 = x(k,:)';    
        net2 = w1k1*O1;
        O2 = 1./(1+exp(-net2));    
        ym(k) = w2k1*O2; 
    
        e(k) = y(k) - ym(k); % 模型误差
        es(num) = es(num) + e(k)^2/2; % 累计误差平方   

        % 训练BP网络
        dw2 = eta*e(k)*O2'; 
        w2 = w2k1 + dw2 + alpha*(w2k1-w2k2); % w2(k)

        df = exp(-net2)./(1+exp(-net2)).^2; % 激励函数的导数
        dw1 = eta*e(k)*w2k1'.*df*O1'; % 矩阵形式计算
        w1 = w1k1 + dw1 + alpha*(w1k1-w1k2); % w1(k)
    
        % 更新数据
        w1k2 = w1k1; w1k1 = w1;
        w2k2 = w2k1; w2k1 = w2;
    end
    eg = es(num);
    if eg <= eg1
        eg1 = eg1 - E/500;
        % fprintf('eg = %f      num = %d\n', eg, num);
    end
end 
toc
figure(1)
plot(t,y,'b',t,ym,'r--')
xlabel('时间t（秒）');
ylabel('实际输出/网络输出');
legend('实际输出', '网络输出','Location','southwest');
figure(2)
plot(1:num,es,'b')
xlabel('训练步数（步）');
ylabel('全局误差 E=0.02');