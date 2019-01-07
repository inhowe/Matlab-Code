% 第二题
clc,clear,close all
p1 = conv([1,0,1],conv([1,3],[1,1]));
p2 = [1,2,1];
[q,r] = deconv(p1,p2)
disp(['商多项式为：',poly2str(q,'t')])
disp(['余多项式为：',poly2str(r,'t')])


