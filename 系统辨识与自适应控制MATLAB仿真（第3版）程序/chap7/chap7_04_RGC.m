% 关系度聚类
clear all; close all;

n = 1; % 聚类算法中的n
N = 15; % 数据向量的个数
rbar = 0.7; % 聚类预设参数：越大分类越细（多）

% 输入输出数据向量
x = [0.6  0.7
    0.1  0.4
    0.5  0.6
    0.3  0.4
    0.7  0.2
    0.6  0.6
    0.9  0.2
    0.7  0.6
    0.2  0.4
    0.6  0.5
    0.8  0.3
    0.2  0.5
    0.8  0.2
    0.2  0.3
    0.8  0.1];

% 计算x与x之间的最大距离
v = x;
for k = 1:N
    for j = k:N
        d(k,j) = norm(x(k,:)-x(j,:));
    end
end
dmax = max(max(d));

% 聚类迭代：Step2-Step5
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

figure(1)
plot(x(:,1),x(:,end),'b*','LineWidth',1); hold on;
plot(c(:,1),c(:,end),'rO','LineWidth',2,'MarkerSize',10); hold off
axis([0 1 0 0.8]);
legend('原始数据','聚类中心','Location','northwest')
xlabel('原始数据x_1');
ylabel('原始数据x_2');