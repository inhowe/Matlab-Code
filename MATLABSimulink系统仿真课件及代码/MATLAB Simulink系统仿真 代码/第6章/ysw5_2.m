clc,clear,close all
ks=[0.4 0.6 0.8];
om=10;
for i=1:3
  num=om*om;
  den=[1 2*ks(i)*om om*om];
  nyquist(num,den);
  axis('square');
  hold on;
  grid on
end

%% 
clc,clear,close all
ks=0.04:0.01:0.707;  
for i=1:length(ks)
  Mr(i)=1/(2*ks(i)*sqrt(1-ks(i)*ks(i)));
end
plot(ks,Mr,'b-');grid;
xlabel('×èÄá±È'),ylabel('Mr');
