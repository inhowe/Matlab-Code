clc,clear,close all
x=1:10;
y=rand(10,1);
plot(x,y,'bo-')
xlabel('x')
ylabel('y')

%%
clc,clear,close all
x=-11:0.1:10;
y=tan(x*pi);
plot(x,y,'r--')
xlabel('x')
ylabel('y')
axis tight

%%
clc,clear,close all
X = -10:1:10; 
Y = -10:1:10;
[X,Y] = meshgrid(X,Y);
Z = - X.^2 - Y.^2 + 10;
surf(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
axis tight
colormap(jet)
shading interp
set(gca,'Ydir','reverse');
set (gcf, 'color', 'w')




