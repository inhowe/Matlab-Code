% 补偿模糊神经网络辨识
clear all; close all;

[t,x] = ode45(@nonsys,[0 20],[1.0 1.0]); % 解微分方程
y = sin(x(:,1)+x(:,2));
[N,n] = size(x); % N为数据个数，n为系统输入的维数

eta = 0.5; % 学习速率
E = 0.02; % 全局误差精度

% 定义输入输出隶属函数中心和宽度及补偿度的初值
M = 25; % 模糊规则数
for j = 1:M
    ak(j,1) = 0.5 + 0.5*floor(j/5);
    ak(j,2) = 0.5 + 0.25*mod(j,5); % 输入隶属函数中心
    for i = 1:n
        sigmak(j,i) = 0.5; % 输入隶属函数宽度
    end
    deltak(j) = 0.5; % 输出隶属函数宽度
    fk(j) = 0.1; hk(j) = 0.1;
    gammak(j) = fk(j)^2/(fk(j)^2+hk(j)^2); % 补偿度
end
% 计算输出隶属函数的中心
x1span = [0.5,0.75,1.25,1.75,2.25,2.5]; 
x2span = [0.4,0.625,0.875,1.125,1.375,1.5]; 
for i1 = 1:5
    for i2 = 1:5
        j = (i1-1)*5 + i2;
        nj = 0;
        for k = 1:N
            if x(k,1)>=x1span(i1) & x(k,1)<= x1span(i1+1)
                if x(k,2)>=x2span(i2) & x(k,2)<= x2span(i2+1)
                    nj = nj + 1;
                    Mj(j,nj) = k;
                end
            end
        end
        bk(j) = 0;
        if nj == 0
            bk(j) = 0.5;
        else
            for i = 1:nj
                bk(j) = bk(j) + y(Mj(j,i));
            end
            bk(j) = bk(j)/nj;
        end
    end
end

eg = 10; % 初始化全局误差
num = 0; % 初始化训练步数
mstep = 100; % 最大训练次数
tic
while (eg >= E) % 由全局误差控制训练次数
% while(num < mstep) % 直接设定训练次数
    num = num+1;
    es(num) = 0;
    for k = 1:N
        % 计算网络输出
        O1 = x(k,1:n);
        bdo = 0;do = 0;
        for j = 1:M
            O3(j) = 1;
            for i = 1:n
                O2(j,i) = exp(-((O1(i)-ak(j,i))/sigmak(j,i))^2);
                O3(j) = O3(j)*O2(j,i);
            end
            gn(j) = 1 - gammak(j) + gammak(j)/n;
            O4(j) = O3(j)^gn(j);
            bdo = bdo + bk(j)*deltak(j)*O4(j);
            do = do + deltak(j)*O4(j); % 公式中的C
        end
        ym(k) = bdo/do;
        e(k) = y(k) - ym(k); 
        es(num) = es(num) + e(k)^2/2; % 累计误差平方  
        
        % 训练网络
        for j = 1:M
            % 对输出隶属函数中心和宽度的训练
            Delta5 = ym(k) - y(k);
            db(j) = Delta5*deltak(j)*O4(j)/do;
            ddelta(j) = Delta5*(bk(j)-ym(k))*O4(j)/do;
            b(j) = bk(j) - eta*db(j);
            delta(j) = deltak(j) - eta*ddelta(j);
            % 对补偿度的训练
            Delta4(j) = Delta5*(bk(j)-ym(k))*deltak(j)/do;
            dgamma(j) = (1/n-1)*Delta4(j)*O4(j)*log(O3(j));
            f(j) = fk(j) - eta*(2*fk(j)*hk(j)^2/(fk(j)^2+hk(j)^2)^2)*dgamma(j);
            h(j) = hk(j) + eta*(2*hk(j)*fk(j)^2/(fk(j)^2+hk(j)^2)^2)*dgamma(j);
            gamma(j) = f(j)^2/(f(j)^2+h(j)^2);
            % 对输入隶属函数中心和宽度的训练
            for i = 1:n
                da(j,i) = 2*Delta5*(bk(j)-ym(k))*deltak(j)*gn(j)*O4(j)*(O1(i)-ak(j,i))/(sigmak(j,i)^2*do);
                dsigma(j,i) = 2*Delta5*(bk(j)-ym(k))*deltak(j)*gn(j)*O4(j)*(O1(i)-ak(j,i))^2/(sigmak(j,i)^3*do);
                a(j,i) = ak(j,i) - eta*da(j,i);
                sigma(j,i) = sigmak(j,i) - eta*dsigma(j,i);
            end
        end
        
        %更新数据
        ak = a; 
        sigmak = sigma;
        fk = f; hk = h;
        gammak = gamma;
        bk = b;
        deltak = delta;
    end
    eg = es(num);
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
ylabel('全局误差 E=0.001');