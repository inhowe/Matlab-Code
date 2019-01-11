clc,clear,close all
g=tf(1,[1 1]); 
nichols(g);
grid on
axis([-135,0,-40,10])

%%
clc,clear,close all
g=tf([1,0],[1]);
nichols(g);
nyquist(g);
grid on

%% 积分环节
clc,clear,close all
g=tf([0,1],[1,0]);
nichols(g);
nyquist(g);
grid on

%%
clc,clear,close all
g=tf(1,[1 1]);
nyquist(g);
% nichols(g);
hold on
g=tf(1,[1 -1]);
nyquist(g,'r');
axis('square');
grid;  

%% 一阶复合微分环节
clc,clear,close all
g=tf([10,1],[0 1]);
nyquist(g);
% nichols(g);
grid on;  
hold on
g=tf([10,-1],[0 1]);
nyquist(g);
% nichols(g);
axis('square');

















