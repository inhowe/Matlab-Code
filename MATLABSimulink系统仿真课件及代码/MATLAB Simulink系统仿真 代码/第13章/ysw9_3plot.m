close all;
clc,
figure(1);
plot(t,y(:,1),'k',t,y(:,2),'r:','linewidth',2);
xlabel('time(s)');ylabel('位置跟踪');
legend('实际信号','仿真结果');
figure(2);
plot(t,cos(t),'k',t,y(:,3),'r:','linewidth',2);
xlabel('time(s)');ylabel('速度追踪');
legend('实际信号','仿真结果');
figure(3);
plot(t,ut(:,1),'r','linewidth',2);
xlabel('time(s)');ylabel('控制输入');