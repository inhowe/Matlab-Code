clc,clear,close all
t=0:0.1:1.5;
Vx=2*t;
Vy=2*t.^2;
Vz=6*t.^3-t.^2;
x=t.^2;
y=(2/3)*t.^3;
z=(6/4)*t.^4-(1/3)*t.^3;  %由速度得到曲线
plot3(x,y,z,'r.-'),       %画飞行轨迹
hold on             
%算数值梯度，也就是重新计算数值速度矢量，这只是为了编程的方便，不是必须的
Vx=gradient(x);
Vy=gradient(y);
Vz=gradient(z);
quiver3(x,y,z,Vx,Vy,Vz),  %画速度矢量图
grid on       % 栅格化
xlabel('x')
ylabel('y')
zlabel('z')

