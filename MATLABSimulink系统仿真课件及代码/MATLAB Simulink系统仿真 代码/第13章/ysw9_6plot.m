close all;
clc
figure(1);
plot(t,y(:,1),'k',t,y(:,2),'r:','linewidth',2);
xlabel('time(s)');ylabel('位置跟踪');
legend('实际信号','仿真结果');
figure(2);
plot(t,y(:,1)-y(:,2),'k','linewidth',2);
xlabel('time(s)');ylabel('位置跟踪误差');
legend('跟踪误差');
figure(3);
plot(t,ut(:,1),'k','linewidth',2);
xlabel('time(s)');ylabel('控制输入');