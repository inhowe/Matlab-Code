clc,clear,close all
X=[1 2 3 4 5];
T=[1 3 5 7 9];
net=newff(X,T);
net = train(net,X,T);
TT = sim(net,X)
gensim(net,-1)


