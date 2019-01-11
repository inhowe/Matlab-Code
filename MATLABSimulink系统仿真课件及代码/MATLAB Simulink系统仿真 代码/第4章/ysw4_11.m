clc,clear,close all
bdclose all
open_system('ysw4_10') % 打开simulink库窗口
A=gcb
B=gcb('ysw4_10')

%%
clc,clear,close all
gcbh

%%
clc,clear,close all
gcs

%%
clc,clear,close all
bdclose all
open_system('ysw4_10') % 打开simulink库窗口
replace_block('built-in/ysw4_10','Scope','Integrator')
replace_block('ysw4_10','Scope','Integrator')

%%
clc,clear,close all
bdclose all
open_system('ysw4_10') % 打开simulink库窗口
get_param('ysw4_10/Scope1','Ymin')
blks = find_system(gcs,'Type','block');
get_param(blks,'BlockType')
get_param('ysw4_10/Sine Wave','DialogParameters')

%%
clc,clear,close all
bdclose all
open_system('ysw4_10') % 打开simulink库窗口
set_param('ysw4_10','Solver','ode15s','StopTime','3000')
set_param('ysw4_10/Sine Wave','Sample time','0.01')

set_param('ysw4_10/Sine Wave','Position',[120,100,150,130])






