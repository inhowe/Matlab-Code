clc,clear,close all
open_system('ysw4_1.slx');
get_param('ysw4_1/Sine Wave','Sample time')

%%
clc,clear,close all
% open_system('ysw4_1.slx');
ysw1 = sprintf('\n');
get_param(['ysw4_1.slx/Signal',ysw1,'Generator'],'Amplitude')

get_param('ysw4_1.slx/Integrator//Noise','Location')