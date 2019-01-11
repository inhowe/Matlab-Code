clc,close all
figure(1);
l = length(simout1);
t = 0:10/(l-1):10;
plot(t,simout1(:,1),'r','linewidth',2)
hold on
plot(t,simout1(:,2),'g','linewidth',2)
plot(t,simout1(:,3),'b','linewidth',2)
legend('定前后轮比例控制的4WS系统','2WS系统','横摆角速度反馈的4WS系统')
figure(2),
plot(t,simout2(:,1),'r','linewidth',2)
hold on
plot(t,simout2(:,2),'b','linewidth',2)
legend('横摆角速度反馈的4WS系统','2WS系统')
%%
clc,close all
figure(1);
l = length(simout1);
t = 0:10/(l-1):10;
plot(t,simout1(:,1),'r','linewidth',2)
hold on
plot(t,simout1(:,2),'g','linewidth',2)
plot(t,simout1(:,3),'b','linewidth',2)
legend('定前后轮比例控制的4WS系统','2WS系统','横摆角速度反馈的4WS系统')
figure(2),
plot(t,simout2(:,1),'r','linewidth',2)
hold on
plot(t,simout2(:,2),'b','linewidth',2)
legend('横摆角速度反馈的4WS系统','2WS系统')




