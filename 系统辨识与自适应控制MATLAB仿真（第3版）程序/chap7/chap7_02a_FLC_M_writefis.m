%非线性系统模糊逻辑控制（建立FIS + writefis）
clear all; close all;

a = newfis('fuzcon_m'); %新建FIS

%添加输入变量e、ec、uc及其隶属函数（三角形）
a = addvar(a, 'input', 'e', [-6 6]);
a = addmf(a, 'input', 1, 'NB', 'trimf', [-8 -6 -4]);
a = addmf(a, 'input', 1, 'NM', 'trimf', [-6 -4 -2]);
a = addmf(a, 'input', 1, 'NS', 'trimf', [-4 -2  0]);
a = addmf(a, 'input', 1, 'Z', 'trimf', [-2  0  2]);
a = addmf(a, 'input', 1, 'PS', 'trimf', [0   2  4]);
a = addmf(a, 'input', 1, 'PM', 'trimf', [2   4  6]);
a = addmf(a, 'input', 1, 'PB', 'trimf', [4   6  8]);

a = addvar(a, 'input', 'ec', [-6 6]);
a = addmf(a, 'input', 2, 'NB', 'trimf', [-8 -6 -4]);
a = addmf(a, 'input', 2, 'NM', 'trimf', [-6 -4 -2]);
a = addmf(a, 'input', 2, 'NS', 'trimf', [-4 -2  0]);
a = addmf(a, 'input', 2, 'Z', 'trimf', [-2  0  2]);
a = addmf(a, 'input', 2, 'PS', 'trimf', [0   2  4]);
a = addmf(a, 'input', 2, 'PM', 'trimf', [2   4  6]);
a = addmf(a, 'input', 2, 'PB', 'trimf', [4   6  8]);

a = addvar(a, 'output', 'uc', [-6 6]);
a = addmf(a, 'output', 1, 'NB', 'trimf', [-8 -6 -4]);
a = addmf(a, 'output', 1, 'NM', 'trimf', [-6 -4 -2]);
a = addmf(a, 'output', 1, 'NS', 'trimf', [-4 -2  0]);
a = addmf(a, 'output', 1, 'Z', 'trimf', [-2  0  2]);
a = addmf(a, 'output', 1, 'PS', 'trimf', [0   2  4]);
a = addmf(a, 'output', 1, 'PM', 'trimf', [2   4  6]);
a = addmf(a, 'output', 1, 'PB', 'trimf', [4   6  8]);

figure(1)
plotmf(a, 'input', 1); %绘制指定变量的所欲隶属函数
xlabel('e、ec、uc');

%建立规则库
r0=[1 1 2 2 3 3 4
    1 2 2 3 3 4 5
    2 2 3 3 4 5 5
    2 3 3 4 5 5 6
    3 3 4 5 5 6 6
    3 4 5 5 6 6 7
    4 5 5 6 6 7 7]; %控制规则表
for i = 1: 7
    for j = 1:7
        r1((i-1)*7+j,:) = [i, j, r0(i,j)];
    end
end
rulelist = [r1, ones(49,2)];
a = addrule(a, rulelist); %添加规则

writefis(a,'fuzcon_m'); %将建好的FIS保存到磁盘的当前目录中