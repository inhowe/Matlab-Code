clc,clear,close all
x=0:0.1:10;
figure('color',[1,1,1])
y=gaussmf(x,[2 5]);
plot(x,y,'m','linewidth',2)
hold on
y1=gaussmf(x,[1 5]);
plot(x,y1,'r','linewidth',2)
y2=gaussmf(x,[1 3]);
plot(x,y2,'b','linewidth',2)
y3=gaussmf(x,[-1 2]);
plot(x,y3,'g','linewidth',2)
y4=gaussmf(x,[5 5]);
plot(x,y4,'k','linewidth',2)
xlabel('gaussmf')



