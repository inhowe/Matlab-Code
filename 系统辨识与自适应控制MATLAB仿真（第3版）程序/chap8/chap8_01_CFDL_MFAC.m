% 基于CFDL的无模型自适应控制
clear all; close all;

ny = 3; nu = 2; % 系统结构参数

N = 1000; % 仿真长度
uk = zeros(nu,1); % 控制输入初值：uk(i)表示u(k-i);
yk = zeros(ny,1); % 系统输出初值
duk = 0; % 控制输入增量初值
yr(1) = 0.5*(-1)^round(0/100); % 期望输出初值

% 设置控制器参数
phihk = 0.5; phih0 = phihk;
eta = 1;
mu = 1;
rho = 0.6;
lambda = 1;
epsilon = 10^(-5);

for k = 1:N
    time(k) = k;
    
    % 系统输出
    if k <= N/2
        y(k) = yk(1)/(1+yk(1)^2) + uk(1)^3;
    else
        y(k) = (yk(1)*yk(2)*yk(3)*uk(2)*(yk(3)-1)+uk(1))/(1+yk(2)^2+yk(3)^2);
    end
    
    % 期望输出
    if k <= N/3
        yr(k+1) = 0.5*(-1)^round(10*k/N);
    elseif k>N/3 & k<=2*N/3
        yr(k+1) = 0.5*sin(10*k*pi/N) + 0.3*cos(20*k*pi/N);
    else
        yr(k+1) = 0.5*(-1)^round(10*k/N);
    end
    
    % 参数估计
    dy(k) = y(k) - yk(1);
    phih(k) = phihk + eta*duk*(dy(k)-phihk*duk)/(mu+duk^2);
    if abs(phih(k))<=epsilon | abs(duk)<=epsilon | sign(phih(k))~=sign(phih0)
        phih(k) = phih0;
    end
    
    % 控制量计算
    u(k) = uk(1) + rho*phih(k)*(yr(k+1)-y(k))/(lambda+phih(k)^2);
    
    % 更新数据
    phihk = phih(k);
    
    duk = u(k) - uk(1);
    for i = nu:-1:2
        uk(i) = uk(i-1);
    end
    uk(1) = u(k);
    
    for i = ny:-1:2
        yk(i) = yk(i-1);
    end
    yk(1) = y(k);
end
figure(1);
plot(time,yr(1:N),'r--',time,y,'b');
xlabel('k'); ylabel('输出跟踪性能');
legend('y_r(k)','y(k)','Location','southeast'); %axis([0 N -1 1]);
figure(2)
plot(time,u,'b');
xlabel('k'); ylabel('控制输入'); %axis([0 N -1 1]);
figure(3)
plot(time,phih,'b');
xlabel('k'); ylabel('伪偏导数估计值'); %axis([0 N 0.4 0.6]);