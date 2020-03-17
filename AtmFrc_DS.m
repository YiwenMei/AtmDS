function [Tad,Tdd,Pad,qd,RHd]=AtmFrc_DS(Ta,LR,Td,LRd,Pa,Z,Zd)
%% Check the inputs
narginchk(7,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Ta',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Ta'));
addRequired(ips,'LR',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'LR'));
addRequired(ips,'Td',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Td'));
addRequired(ips,'LRd',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'LRd'));
addRequired(ips,'Pa',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Pa'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Z'));
addRequired(ips,'Zd',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Zd'));
parse(ips,Ta,LR,Td,LRd,Pa,Z,Zd);
clear ips

%% Downscaling
Z=readCls(Z);
Zd=readCls(Zd);
dZ=Zd-imresize(Z,size(Zd),'bilinear');
clear Z Zd

% Air temperature (K)
LR=readCls(LR);
Ta=readCls(Ta);
Tad=imresize(Ta,size(dZ),'bilinear')+imresize(LR,size(dZ),'bilinear').*dZ;

% Air pressure (Pa)
Pa=readCls(Pa);
Pad=Pair_Adj(imresize(Pa,size(dZ),'bilinear'),0,imresize(Ta,size(dZ),'bilinear'),dZ,LR);
clear Pa LR

% Dew point temperature (K)
LRd=readCls(LRd);
Td=readCls(Td);
Td(Td>Ta)=Ta(Td>Ta); % Set Td > Ta to Ta
Tdd=imresize(Td,size(dZ))+imresize(LRd,size(dZ),'bilinear').*dZ;
clear Td LRd Ta
Tdd(Tdd>Tad)=Tad(Tdd>Tad); % Set Td > Ta to Ta

% Humidity
[RHd,qd]=Hum_Cal(Tad,Pad,'DewPoint',Tdd);
end

function v2d=readCls(vb)
if isa(vb,'char')
  v2d=matfile(vb);
  vb=cell2mat(who(v2d));
  eval(sprintf('v2d=v2d.%s;',vb));
else
  v2d=vb;
end
end
