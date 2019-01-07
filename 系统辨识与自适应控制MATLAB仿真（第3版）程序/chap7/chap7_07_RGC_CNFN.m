% 基于关系度聚类的补偿模糊神经网络辨识
clear all; close all;

[t,x] = ode45(@nonsys,[0 20],[1.0 1.0]); % 解微分方程
y = sin(x(:,1)+x(:,2));
[N,n] = size(x); % N为数据个数，n为系统输入的维数

eta = 0.5; % 学习速率
E = 0.02; % 全局误差精度

rbar = 0.932; % 聚类预设参数：越大分类越细（多）

% 计算x与x之间的最大距离
x = [x, y]; % 输入输出数据
v = x;
for k = 1:N
    for j = k:N
        d(k,j) = norm(x(k,:)-x(j,:));
    end
end
dmax = max(max(d));

% 聚类迭代：聚类算法的Step2-Step5
w = zeros(size(v));
cstep = 0; % 聚类迭代步数
while (1)
    for k = 1:N
        rvsum = 0;
        rsum = 0;
        for j = 1:N
            r(k,j) = 1-norm(v(k,:)-v(j,:))/dmax; 
            if r(k,j) < rbar
                r(k,j) = 0;
            end
            rvsum = rvsum + r(k,j)*v(j,:);
            rsum = rsum + r(k,j);
        end
        w(k,:) = rvsum/rsum;
    end
    cstep = cstep + 1;
    wvsum = sum(sum(roundn(w,-10)==roundn(v,-10)));
    if wvsum == N*(n+1)
        break;
    else
        v = w;
    end
end

% 从v中取出聚类中心
M = 1; % M为聚类数
c(1,:) = v(1,:); % c为聚类的中心
for k = 2:N
    for j = 1:M
        vsum(j) = sum(roundn(v(k,:),-10)==roundn(c(j,:),-10));
    end
    if max(vsum) ~= n+1 
        M = M + 1;
        c(M,:) = v(k,:);
    end
end
ak = c(:,1:n); % 输入隶属函数中心的初值
bk = c(:,n+1); % 输出隶属函数中心的初值

% 根据聚类结果，将原始数据归类
nn = zeros(M,1);
for k = 1:N
    for j = 1:M
        if sum(abs(v(k,:)-c(j,:))) < 10^(-10)
            nn(j) = nn(j) + 1;
            xx(nn(j),:,j) = x(k,:);
            break;
        end
    end
end

% 设置输入输出隶属函数宽度的初值
mx(1,:) = 0.1*ones(1,n+1);
for j = 1:M
    if nn(j) ~= 1
        mx(j,:) = max(abs(xx(1:nn(j),:,j)-ones(nn(j),1)*c(j,:)));
    else
        mx(j,:) = mean(mx);
    end
    sigmak(j,:) = 2*mx(j,1:n);
    deltak(j,1) = 2*mx(j,end);
end
% 调整非常小的隶属函数宽度
for j = 1:M
    for i = 1:n
        if sigmak(j,i) < 0.1
            sigmak(j,i) = mean(sigmak(:,i));
        end
    end
    if deltak(j) < 0.1
        deltak(j) = mean(deltak);
    end
end

% 设置补偿度的初值
for j = 1:M
    fk(j) = 0.1; hk(j) = 0.1;
    gammak(j) = fk(j)^2/(fk(j)^2+hk(j)^2); % 补偿度
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
        O1 = x(k,:);
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
plot(t,y,'k',t,ym,'k--')
xlabel('时间t（秒）');
ylabel('实际输出/网络输出');
legend('实际输出', '网络输出','Location','southwest');
figure(2)
plot(1:num,es,'k')
xlabel('训练步数（步）');
ylabel('全局误差 E=0.02');