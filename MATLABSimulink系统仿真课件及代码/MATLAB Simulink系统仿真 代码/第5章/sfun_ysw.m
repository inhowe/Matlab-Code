function [sys,x0,str,ts,simStateCompliance] = sfun_ysw(t,x,u,flag)
switch flag,
  % Initialization %
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  % Derivatives %
  case 1,
    sys=mdlDerivatives(t,x,u);
  % Update %
  case 2,
    sys=mdlUpdate(t,x,u);
  % Outputs %
  case 3,
    sys=mdlOutputs(t,x,u);
  % GetTimeOfNextVarHit %
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);
  % Terminate %
  case 9,
    sys=mdlTerminate(t,x,u);
  % Unexpected flags %
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

% mdlInitializeSizes
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;   % at least one sample time is needed
sys = simsizes(sizes);
% initialize the initial conditions
x0  = [];
% str is always an empty matrix
str = [];
% initialize the array of sample times
ts  = [];
simStateCompliance = 'UnknownSimState';

% mdlDerivatives
function sys=mdlDerivatives(t,x,u)
sys = [];

% mdlUpdate
function sys=mdlUpdate(t,x,u)
sys = [];

% mdlOutputs
function sys=mdlOutputs(t,x,u)
sys = 2*u;

% mdlGetTimeOfNextVarHit
function sys=mdlGetTimeOfNextVarHit(t,x,u)
sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% mdlTerminate
function sys=mdlTerminate(t,x,u)
sys = [];
