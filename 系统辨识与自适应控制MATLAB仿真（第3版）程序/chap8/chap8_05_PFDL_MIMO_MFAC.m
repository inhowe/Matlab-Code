% 基于PFDL的MIMO无模型自适应控制
clear all; close all;

m = 2; % m输入m输出系统
ny = 2; nu = 2; % 系统结构参数

N = 800; % 仿真长度
L = 2; % 控制输入线性化长度常数
uk = zeros(m,nu); % 控制输入初值：uk(1,i)表示u1(k-i);
yk = zeros(m,ny); % 系统输出初值
duk = zeros(m,L); % 控制输入增量初值
yr(1,1) = 5*sin(0/50) + 2*cos(0/20); % 期望输出初值
yr(2,1) = 2*sin(0/50) + 5*cos(0/20);

% 设置控制器参数
Phihk = [[1.5  0.1; 0.1  1.5], zeros(m,m*(L-1))]; Phih0 = Phihk(:,1:m);
eta = 0.5;
mu = 1;
rho = 0.5*ones(L,1);
lambda = 0.01;
b1 = 0.3;
b2 = 1.3;
alpha = 1.5;

for k = 1:N
    time(k) = k;
    
    % 系统输出
    y(1,k) =  -0.8*sin(yk(1,1)) - 0.16*yk(1,2)...
             + uk(1,1) + 1.7*uk(1,2) - 0.5*uk(2,2) + uk(2,2)/(1+uk(2,2)^2);
    y(2,k) = -0.8*yk(2,1) - 0.16*yk(2,2)...
             + uk(2,1) + 2*uk(2,2) + 0.3*uk(1,2) + uk(1,2)/(1+uk(1,2)^2);

    % 期望输出
    yr(1,k+1) = 2*sin(k/50) + 5*cos(k/20);
    yr(2,k+1) = 5*sin(k/50) + 2*cos(k/20);
    
    % 参数估计
    dy(:,k) = y(:,k) - yk(:,1); 
    dHk = [];
    for i = 1:L
        dHk = [dHk; duk(:,i)];
    end
    Phih(:,:,k) = Phihk + eta*(dy(:,k)-Phihk*dHk)*dHk'/(mu+norm(dHk)^2);
    
    for i = 1:m
        for j = 1:m
            if i==j
                if abs(Phih(i,j,k))<b2 | abs(Phih(i,j,k))>alpha*b2 | sign(Phih(i,j,k))~=sign(Phih0(i,j))
                    Phih(i,j,k) = Phih0(i,j);
                end
            else
                if abs(Phih(i,j,k))>b1 | sign(Phih(i,j,k))~=sign(Phih0(i,j))
                    Phih(i,j,k) = Phih0(i,j);
                end
            end
        end
    end
    for p = 1:L
        for i = 1:m
            for j = 1:m
                phih((p-1)*m*m+(i-1)*m+j,k) = Phih(i,(p-1)*m+j,k);
            end
        end
    end
    
    % 控制量计算
    sumrpdu = zeros(m,1);
    for i = 2:L
        sumrpdu = sumrpdu + rho(i)*Phih(:,(i-1)*m+1:i*m,k)*duk(:,i-1);
    end
    du(:,k) = Phih(:,1:m,k)'*(rho(1)*(yr(:,k+1)-y(:,k)) - sumrpdu)/(lambda+norm(Phih(:,1:m,k))^2);
    u(:,k) = uk(:,1) + du(:,k);
    
    % 更新数据
    Phihk = Phih(:,:,k);
    
    for i = L:-1:2
        duk(:,i) = duk(:,i-1);
    end
    duk(:,1) = du(:,k);
    
    for i = nu:-1:2
        uk(:,i) = uk(:,i-1);
    end
    uk(:,1) = u(:,k);
    
    for i = ny:-1:2
        yk(:,i) = yk(:,i-1);
    end
    yk(:,1) = y(:,k);
end
figure(1);
plot(time,yr(1,1:N),'r--',time,y(1,:),'b');
xlabel('k'); ylabel('输出跟踪性能'); %axis([0 N -8 8]);
legend('y_{r1}(k)','y_1(k)','Location','best');
figure(2);
plot(time,yr(2,1:N),'r--',time,y(2,:),'b');
xlabel('k'); ylabel('输出跟踪性能'); %axis([0 N -8 8]);
legend('y_{r2}(k)','y_2(k)','Location','best');
figure(3)
plot(time,u);
xlabel('k'); ylabel('控制输入'); %axis([0 N -5 4]);
legend('u_1(k)','u_2(k)','Location','best');
if L==1
    figure(4)
    plot(time,phih(1:m*m,:),'LineWidth',1.5);
    % hold on; plot(time,0.3*ones(1,N),time,1.3*ones(1,N),time,1.95*ones(1,N));
    xlabel('k'); ylabel('伪雅克比矩阵估计值'); %axis([0 N 0 2]);
    legend('\phi_{11}(k)估计值','\phi_{12}(k)估计值','\phi_{21}(k)估计值','\phi_{22}(k)估计值','Location','best');
elseif L==2
    figure(4)
    plot(time,phih(1:m*m,:),'LineWidth',1.5);
    % hold on; plot(time,0.3*ones(1,N),time,1.3*ones(1,N),time,1.95*ones(1,N));
    xlabel('k'); ylabel('伪雅克比矩阵估计值'); %axis([0 N 0 2]);
    legend('\phi_{111}(k)估计值','\phi_{121}(k)估计值','\phi_{211}(k)估计值','\phi_{221}(k)估计值','Location','best');
    figure(5)
    plot(time,phih(m*m+1:end,:),'LineWidth',1.5);
    % hold on; plot(time,0.3*ones(1,N),time,1.3*ones(1,N),time,1.95*ones(1,N));
    xlabel('k'); ylabel('伪雅克比矩阵估计值'); %axis([0 N 0 2]);
    legend('\phi_{112}(k)估计值','\phi_{122}(k)估计值','\phi_{212}(k)估计值','\phi_{222}(k)估计值','Location','best');
end