close all;
clc
figure(1);
plot(t,y(:,1),'r',t,y(:,3),'b','linewidth',2);
xlabel('time(s)');ylabel('位置跟踪');
legend('实际信号','仿真结果');
figure(2);
plot(t,y(:,2),'r',t,y(:,4),'b','linewidth',2);
xlabel('time(s)');ylabel('速度跟踪');
legend('实际信号','仿真结果');
figure(3);
plot(t,u(:,1),'r','linewidth',2);
xlabel('time(s)');ylabel('控制输入');