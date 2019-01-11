clc,clear,close all
ww1=0.1:0.01:10;
for i=1:length(ww1)
  Lw=(-20)*log10(sqrt(1+ww1(i)^2));
  if ww1(i)<=1 Lw1=0;
  else Lw1=(-20)*log10(ww1(i));
  end
  m(i)=Lw-Lw1;
end
ab=semilogx(ww1,m,'b-');
set(ab,'LineWidth',2);grid;
xlabel('w/w1'),ylabel('Îó²î/dB');

%%
clc,clear,close all
ks=[0.1 0.2 0.3 0.5 0.7 1.0];
om=10;
for i=1:length(ks)
    num=om*om;
    den=[1 2*ks(i)*om om*om];         
    bode(num,den);
    hold on;      
end
grid;


%%
clc,clear,close all
ks=[0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.8 1.0];
wwn=0.1:0.01:10;
for i=1:length(ks)
     for k=1:length(wwn)
         Lw=-20*log10(sqrt((1-wwn(k)^2)^2+(2*ks(i)*wwn(k))^2));
         if wwn(k)<=1 Lw1=0;
         else Lw1=-40*log10(wwn(k));
         end
         m(k)=Lw-Lw1;
     end
     ab=semilogx(wwn,m,'b-');
     set(ab,'linewidth',1.5);hold on;
end
grid;


%%
clc,clear,close all
num=[64,128];
a1 = conv([1,0], [1,0.5]);
a2 = conv(a1, [1,3.2,64]);
den=[a2];         
bode(num,den);
hold on;      
grid;




