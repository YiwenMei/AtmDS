function To=Temp_Adj(Ti,Hi,varargin)
%% Check the inputs
narginchk(2,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Ti',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Ti'));
addRequired(ips,'Hi',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Hi'));

addOptional(ips,'Ho',0,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Ho'));
addOptional(ips,'LR',-0.0065,@(x) validateattributes(x,{'double','char'},{'nonempty'},...
    mfilename,'LR'));

parse(ips,Ti,Hi,varargin{:});
Ho=ips.Results.Ho;
LR=ips.Results.LR;
clear ips varargin

%% Adjust temperature
Ti=readCls(Ti);
Hi=readCls(Hi);
LR=readCls(LR);
LR=imresize(LR,size(Ti),'bilinear');

To=LR.*(Ho-Hi)+Ti; % Adjust temperature
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
