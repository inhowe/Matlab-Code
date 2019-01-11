%解析法辨识二阶传函

clc
close all
clear all
%实验数据，数据来源： 《系统辨识方法及应用》 .国防工业出版社
% t=[1 3 5 7 9 11 13 15 17 19];
% y=[0.149086 0.5890067 0.830617 0.933990 0.973980 0.991095 0.995868 0.998680 0.999490 0.999850]; 
t=0:0.1:10;
t=t';
num=[24];
den=[1 10 24];
y=step(num,den,t);

y2=log(1-y);
plot(t,y2,'*');
grid on
pm=polyfit(t,y2,1)
value=polyval(pm,t);
hold on
plot(t,value,'r')
title('\fontname{ 黑体 }\fontsize{20}y(t)=at+b')
w2=-pm(1)
w1=w2/(1-exp(-pm(2)))
T=1/sqrt(w1*w2)
theta=(w1+w2)/(2*sqrt(w1*w2))
z=[];
p=[-w1 -w2];
k=w1*w2;
sys=zpk(z,p,k)
figure
step(sys);
% axis([0 20 0 1.2])
hold on
% plot(t,y,'r*')
step(num,den)
zpk(tf(num,den))