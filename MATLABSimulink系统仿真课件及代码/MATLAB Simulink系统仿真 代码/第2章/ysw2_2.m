% PID控制阶跃信号代码

clc,clear,close all
% PID参数 --- 可手动调节
kp = 0.4267;
ki = 7.7329;
kd = 1.607;

iter=100;     % 增量式PID迭代步数
errori = 0;
errord_1 = 0;
u2=1;
y(1)=0;
a=0.05;
for k=1:iter

    error = u2 - y(k); 
    errori = errori+error;
    if errori>=a
        errori=a;
    end
    if errori<=-a
        errori=-a;
    end
    errord = (error - errord_1)*a;
    y(k+1) = kp*error+ki*errori+kd*errord;
    errord_1=error;

end
plot([1,iter],[u2,u2],'r-','linewidth',2)
hold on
plot(1:iter,y(2:end),'linewidth',2)
axis([0,iter,0,1.4])