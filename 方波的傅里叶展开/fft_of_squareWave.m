% 方波的傅里叶变换
% 频率为fHz的方波的符号定义法：squareWave=sign(sin(2*pi*f*t)
% 方波的傅里叶展开公式：squareWave=4/pi*{sin(2*pi*f*t)+1/3*sin(6*pi*f*t)+1/5*sin(10*pi*f*t)+...}
clc;clear;close all;
t=0:0.01:5;
freq=1;
omiga=2*pi*freq;
squareWave=sign(sin(omiga*t));

fft=4/pi*(sin(omiga*t));
subplot(3,1,1);
plot(squareWave);hold on
plot(fft);

fft=4/pi*(sin(omiga*t)+1/3*sin(3*omiga*t));
subplot(3,1,2);
plot(squareWave);hold on
plot(fft);

fft=4/pi*(sin(omiga*t)+1/3*sin(3*omiga*t)+1/5*sin(5*omiga*t));
subplot(3,1,3);
plot(squareWave);hold on
plot(fft);
