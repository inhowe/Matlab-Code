clc,clear,close all
Z1=[1 2 4;3 4 1];
Z2=[1:3;2:4];    
b=[0;1];
q=4;                       
Z=concur(b,q) 
X1=netsum(Z1,Z2),
X2=netprod(Z1,Z2) %计算向量的和与积
