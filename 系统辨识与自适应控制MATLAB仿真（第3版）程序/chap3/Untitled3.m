cla
G=tf([3],[1 2 3]);
Gz=c2d(G,0.1,'i');
Gz=tf([3 2 0 0],[1 -1 -0.1 0.8],0.1);
% step(Gz)
% step(feedback(Gz,1))
am=[1;-conv([1 -0.8],[1 0.4])']'
bm=[3 2 0 0]; %参考模型参数（参考模型中含有yr(k)，注意nb的使用！）
Gz=tf(bm,am,1);
[A B C D]=tf2ss(am,bm)
syms z
H=D+C*inv(z*eye(3)-A)*B
eig(H)
