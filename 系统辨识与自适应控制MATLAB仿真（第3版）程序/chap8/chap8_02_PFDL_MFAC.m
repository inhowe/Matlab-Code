% 基于PFDL的无模型自适应控制
clear all; close all;

ny = 3; nu = 3; % 系统结构参数

N = 800; % 仿真长度
L = 3; % 控制输入线性化长度常数
uk = zeros(nu,1); % 控制输入初值：uk(i)表示u(k-i);
yk = zeros(ny,1); % 系统输出初值
dhk = zeros(L,1); % 控制输入增量初值
yr(1) = 5*sin(0/50) + 2*cos(0/20); % 期望输出初值

% 设置控制器参数
Phihk = [2; zeros(L-1,1)]; Phih0 = Phihk;
eta = 0.5;
mu = 1;
rho = 0.5*ones(L,1);
lambda = 0.01;
epsilon = 10^(-5);

for k = 1:N
    time(k) = k;
    
    % 系统输出
    if k <= 400
        y(k) = 2.5*yk(1)*yk(2)/(1+yk(1)^2+yk(2)^2) + 1.2*uk(1) + 0.09*uk(1)*uk(2) + 1.6*uk(3) + 0.7*sin(0.5*(yk(1)+yk(2)))*cos(0.5*(yk(1)+yk(2)));
    else
        y(k) = 5*yk(1)*yk(2)/(1+yk(1)^2+yk(2)^2+yk(3)^2) + uk(1) + 1.1*uk(2);
    end
    
    % 期望输出
    yr(k+1) = 5*sin(k/50) + 2*cos(k/20);
    
    % 参数估计
    dy(k) = y(k) - yk(1);
    Phih(:,k) = Phihk + eta*dhk*(dy(k)-Phihk'*dhk)/(mu+norm(dhk)^2);
    if norm(Phih(:,k))<=epsilon | norm(dhk)<=epsilon | sign(Phih(1,k))~=sign(Phih0(1))
        Phih(:,k) = Phih0;
    end
    
    % 控制量计算
    sumrpdu = 0;
    for i = 2:L
        sumrpdu = sumrpdu + rho(i)*Phih(i,k)*dhk(i-1);
    end
    du(k) = rho(1)*Phih(1,k)*(yr(k+1)-y(k))/(lambda+Phih(1,k)^2) - Phih(1,k)*sumrpdu/(lambda+Phih(1,k)^2);
    u(k) = uk(1) + du(k);
    
    % 更新数据
    Phihk = Phih(:,k);
    
    for i = L:-1:2
        dhk(i) = dhk(i-1);
    end
    dhk(1) = du(k);
    
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
legend('y_r(k)','y(k)','Location','southeast'); %axis([0 N -8 8]);
figure(2)
plot(time,u,'b');
xlabel('k'); ylabel('控制输入'); %axis([0 N -5 3]);
figure(3)
plot(time,Phih);
xlabel('k'); ylabel('伪梯度估计值'); %axis([0 N -2 3]);
if L==3
    legend('\phi_1(k)估计值','\phi_2(k)估计值','\phi_3(k)估计值');
end