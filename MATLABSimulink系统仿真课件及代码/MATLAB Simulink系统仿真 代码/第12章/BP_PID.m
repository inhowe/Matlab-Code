function yout = BP_PID(u1,rin)

persistent x u_1 u_2 u_3 u_4 u_5 y_1 y_2 y_3
persistent xite alfa IN H Out wi wi_1 wi_2  wi_3 wo wo_1 wo_2 wo_3
persistent Oh error_2 error_1
persistent k
if u1==0
    xite=0.25;
    alfa=0.05;
    IN=4;H=5;Out=3;  %NN Structure
    wi=[-0.2846    0.2193   -0.5097   -1.0668;
        -0.7484   -0.1210   -0.4708    0.0988;
        -0.7176    0.8297   -1.6000    0.2049;
        -0.0858    0.1925   -0.6346    0.0347;
         0.4358    0.2369   -0.4564   -0.1324];
    %wi=0.50*rands(H,IN);
    wi_1=wi;wi_2=wi;wi_3=wi;
    wo=[1.0438    0.5478    0.8682    0.1446    0.1537;
        0.1716    0.5811    1.1214    0.5067    0.7370;
        1.0063    0.7428    1.0534    0.7824    0.6494];
    %wo=0.50*rands(Out,H);
    wo_1=wo;wo_2=wo;wo_3=wo;
    x=[0,0,0];
    u_1=0;u_2=0;u_3=0;u_4=0;u_5=0;
    y_1=0;y_2=0;y_3=0;
    Oh=zeros(H,1);    %Output from NN middle layer
    I=Oh;             %Input to NN middle layer
    error_2=0;
    error_1=0;
    k=2;
end

%Unlinear model
a = 1.2*(1-0.8*exp(-0.1*k));
yout = a*y_1/(1+y_1^2)+u_1;
error = rin-yout;
xi=[rin,yout,error,1];

x(1)=error-error_1;
x(2)=error;
x(3)=error-2*error_1+error_2;
 
epid=[x(1);x(2);x(3)];
I=xi*wi';
for j=1:1:H
    Oh(j)=(exp(I(j))-exp(-I(j)))/(exp(I(j))+exp(-I(j))); %Middle Layer
end
K=wo*Oh;             %Output Layer
for l=1:1:Out
    K(l)=exp(K(l))/(exp(K(l))+exp(-K(l)));        %Getting kp,ki,kd
end
kp=K(1);ki=K(2);kd=K(3);
Kpid=[kp,ki,kd];
 
du=Kpid*epid;
u=u_1+du;
if u>=10       % Restricting the output of controller
   u=10;
end
if u<=-10
   u=-10;
end
 
dyu=sign((yout-y_1)/(u-u_1+0.0000001));
 
%Output layer
for j=1:1:Out
    dK(j)=2/(exp(K(j))+exp(-K(j)))^2;
end
for l=1:1:Out
    delta3(l)=error*dyu*epid(l)*dK(l);
end
 
for l=1:1:Out
   for i=1:1:H
       d_wo=xite*delta3(l)*Oh(i)+alfa*(wo_1-wo_2);
   end
end
    wo=wo_1+d_wo+alfa*(wo_1-wo_2);
%Hidden layer
for i=1:1:H
    dO(i)=4/(exp(I(i))+exp(-I(i)))^2;
end
segma=delta3*wo;
for i=1:1:H
   delta2(i)=dO(i)*segma(i);
end

d_wi=xite*delta2'*xi;
wi=wi_1+d_wi+alfa*(wi_1-wi_2);
 %Parameters Update
u_5=u_4;u_4=u_3;u_3=u_2;u_2=u_1;u_1=u;   
y_2=y_1;y_1=yout;
   
wo_3=wo_2;
wo_2=wo_1;
wo_1=wo;
   
wi_3=wi_2;
wi_2=wi_1;
wi_1=wi;
 
error_2=error_1;
error_1=error;





