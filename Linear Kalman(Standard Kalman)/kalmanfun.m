function [xe,pk,p1]=kalmanfun(A,C,Q,R,xe,z,p)
%This function is to calculate the estimation state by Kalman filter.
%This function is to calculate the estimation state by Kalman filter.
%输入参数：A—过程矩阵    C—测量矩阵    Q—过程噪声方差    R—测量噪声方差
%          xe—前一步状态估计值     x(k-1|k-1)   
%           z—当前测量值           z(k)
%           p—前一步状态估计方差   P(k-1|k-1)
%输出参数：xe—当前步状态估计值     x(k|k)   
%          pk—向前一步状态估计方差 P(k|k-1) 
%          p1—当前步状态估计方差   P(k|k)
   xe=A*xe;
   p1=A*p*A'+Q;
   K=p1*C'*inv(C*p1*C'+R);
   xe=xe+K*(z-C*xe);
   pk=(eye(size(p1))-l*C)*p1;