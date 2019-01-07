clc,clear,close all
bdclose all
new_system('ysw4_9'); % 新建一个ysw4_9系统

%%
clc,clear,close all
bdclose all
open_system('ysw4_4') % 打开simulink库窗口

%%
clc,clear,close all
bdclose all
open_system('ysw4_7') % 打开simulink库窗口
save_system('ysw4_7','ysw4_10')

%%
clc,clear,close all
bdclose all
open_system('ysw4_10') % 打开simulink库窗口
add_line('ysw4_10','Sine Wave/1','Scope/1')
delete_line('ysw4_10','Sine Wave/1','Scope/1')



