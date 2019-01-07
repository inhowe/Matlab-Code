clc,clear,close all
bdclose
open_system('ysw4_7.slx');
add_block('built-in/Sine Wave','ysw4_7/Sine Wave');

%%
clc,clear,close all
bdclose
open_system('ysw4_7.slx');
delete_block('ysw4_7/Sine Wave')
%%
clc,clear,close all
bdclose
open_system('ysw4_7.slx');
bdroot('ysw4_7/Scope')

%%
find_system
%%
open_bd_ysw = find_system('Type','block_diagram')
%%
open_bd_ysw1 = find_system('ysw4_7/Subsystem','SearchDepth',1,'blockType','Abs')
%%
open_bd_ysw1 = find_system('ysw4_7/Subsystem','FindAll','on','type','line')









