function my_fun
    format short
    x0 = [1; 1];           % 初始点
    options=optimset('Display','iter');   % 优化
    [x,fval] = fsolve(@myfun,x0,options) 
end
function F = myfun(x)
F = [2*x(1)+3*x(2)-exp(-2*x(1));
      x(1)-x(2)+exp(-sqrt(x(2)))];
end