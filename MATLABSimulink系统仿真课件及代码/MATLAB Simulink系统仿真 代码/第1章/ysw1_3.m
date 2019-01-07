clc,clear,close all
sin_f=@sin;
a=sin_f(pi)
myadd = @(x,y) x*sin(y);
b1=myadd(1,1)

%%
clc,clear,close all
a=1;
if a==1
    b=0
else
    b=1
end

%%
clc,clear,close all
a=11;
switch a
    case 1
    b=0
    otherwise
    b=1
end



