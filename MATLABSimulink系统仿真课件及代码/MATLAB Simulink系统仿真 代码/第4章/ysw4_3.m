clc,clear,close all
bdclose
new_system('ysw4_4'); % 新建一个ysw4_4系统
open_system('ysw4_4') % 打开simulink库窗口
add_block('built-in/Sine Wave','ysw4_4/Sine Wave');
add_block('built-in/Scope','ysw4_4/Scope');
add_block('built-in/Mux','ysw4_4/Mux');
add_block('built-in/Integrator','ysw4_4/Integrator');
%%
add_line('ysw4_4','Sine Wave/1','Mux/1')
%%
add_line('ysw4_4',[90,20;80,85;120,90])
%%
add_line('ysw4_4','Integrator/1','Mux/2')
%%
add_line('ysw4_4','Mux/1','Scope/1')
