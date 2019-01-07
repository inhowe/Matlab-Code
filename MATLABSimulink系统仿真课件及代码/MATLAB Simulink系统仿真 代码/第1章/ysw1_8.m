clc,clear,close all
X=fminsearch('2*x(1)^3+x(1)*x(2)^4-10*x(1)*x(2)',  [0,0])

%%
clc,clear,close all
f=[1,1,1];
A =[1 -1  1;3  2  4;3  2  0];
b = [20; 42; 30];
lb = zeros(3,1);
[x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb)
lambda.ineqlin
lambda.lower
