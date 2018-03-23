clc;
clear;
close all
subplot(4,1,1)
t=0:0.01:9.99;
m=(t-5).*(t-5)+10;%样例数据
v=2;%方差
s3=m+sqrt(v)*randn(1,1000);%带噪声样例数据
plot(t,s3)
ea=[0;0;0];           %当前估计值
M=diag([1 1 1])*1000;     %M为一大的初值
estiamte=[];        %所有估计值
W=1/v;              %权值矩阵，取方差倒数

xlabel('Press Any Key To Continue!')
pause
for i=1:1000   
    h=[t(i)*t(i) t(i) 1]; %测量矩阵(变量值）at^2+bt+c=0; 三角函数如何表示？
    M=inv(inv(M)+h'*h);
    ea=ea+M*h'*W*(s3(i)-h*ea);%新估计值=旧估计值+增益*新息

    estiamte=[estiamte,ea]; %记录估计值的变化

    subplot(4,1,1)
    hold on
    lastP1=plot(t,ea(1)*t.*t+ea(2)*t+ea(3),'r'); %绘制叠加拟合曲线
    xx=ones(round(ea(1)*t(1000)*t(1000)+ea(2)*t(1000)+ea(3)),1)*i*0.01; %绘制绿线进度条
    yy=1:round(ea(1)*t(1000)*t(1000)+ea(2)*t(1000)+ea(3));
    lastPx=plot(xx,yy,'g');
    str=['剩余',num2str(1000-i),'个点'];
    xlabel(str)
    hold off

    subplot(4,1,2)
    lastP2=plot(estiamte(1,:)); %绘制a系数的变化过程
    subplot(4,1,3)
    lastP3=plot(estiamte(2,:)); %绘制b系数的变化过程
    subplot(4,1,4)
    lastP4=plot(estiamte(3,:)); %绘制c系数的变化过程
    pause(0.01);
    if(i<1000)
        delete(lastP1)
        delete(lastPx)
    end
end
subplot(4,1,1)
xlabel('t')
ylabel('拟合曲线')
subplot(4,1,2)
ylabel('系数a')
subplot(4,1,3)
ylabel('系数b')
subplot(4,1,4)
xlabel('points')
ylabel('系数c')