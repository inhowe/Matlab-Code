clf;
% clear all;
EndTime=1;%仿真截止时间
sampleTime=0.01;%采样间隔
t=0:sampleTime:EndTime;
L=EndTime/sampleTime; %仿真长度

G=tf([0.98407],[13.407 1]);%测量系统
Gz=c2d(G,sampleTime,'i')
yk=0;
yk1=0;
yk2=0;
yk3=0;
uk=ones(100000,1);
uk1=0;
uk2=0;
for k=1:L
%     yk(k)=1*yk1-0.7*yk2-0.5*yk3+0*uk(k)+3*uk(k)+2*uk1;
    yk(k)=[thetape_1(1) 0 0]*[yk1 yk2 yk3]'+[thetape_1(2) thetape_1(3) 0]*[uk(k) uk1 uk2]';
    yk3=yk2;
    yk2=yk1;
    yk1=yk(k);
    uk2=uk1;
    uk1=uk(k);
end
plot(yk)
hold on;
% step(tf([3 2 0 0],[1 -1 0.7 0.5],1),1:1:50)
% step(Gz)
