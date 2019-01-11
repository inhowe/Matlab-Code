% % clear; 
% yk1=0;yk2=0;yk3=0;
% yk=0;yr1=0;yr2=0;
% for t=1:L
% %     yk(t,1)=thetape_1(1:3)'*[yk1 yk2 yk3]'+thetape_1(4:5)'*[yr(t,1) yr1]';
%     yk(t,1)=thetape_1(1:1)'*[yk1]'+thetape_1(2:2)'*[yr(t,1)]';
%     yk3=yk2;
%     yk2=yk1;
%     yk1=yk(t,1);
%     yr1=yr(t,1);
% end
% figure
% plot(yk),hold on,plot(ym)

% clear;
cla
% G=tf([1],[1 1])
G=tf([3],[1 2 3])
Gz=c2d(G,0.1,'i')
seqGz=step(Gz);
yk1=0;yk2=0;yk3=0;
yk=0;yr1=0;yr2=0;
yr=ones(1000,1);
% yr(1)=0;
for t=1:length(seqGz)
    yk(t,1)=0.4*yk1+0.32*yk2+3*yr(t,1)+2;
    yk2=yk1;
    yk1=yk(t,1);
end
% figure
hold on
% yk=[0;yk];
plot(yk,'r')
plot(seqGz)
% ,hold on,plot(ym)