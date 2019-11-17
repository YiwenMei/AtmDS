% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/12/2019

%% Functionality
% This function adjusts air pressure to a user-specified height.

%% Input
% Pi: spatial map class (V2DCls.m) object or workspace variable for air pressure
%      at input height (Pa);
% Hi: V2DCls.m object or workspace variable for input height of air pressure
%      (m above ground);
% Ti: spatial map class (V2DCls.m) object or workspace variable for air temperature
%      at input height (K).

% Hn: V2DCls.m object or workspace variable for output height of air pressure
%      (m a.g.);
% LR: V2DCls.m object or workspace variable for air temperature lapse rate (K/m).

%% Output
% Po: air pressure at output height (Pa).

%% Additional note
% Require V2DCls.m.

function Po=Pair_Adj(Pi,Hi,Ti,varargin)
%% Check the inputs
narginchk(3,5);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Pi',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Pi'));
addRequired(ips,'Hi',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hi'));
addRequired(ips,'Ti',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Ti'));

addOptional(ips,'Ho',0,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Ho'));
addOptional(ips,'LR',-0.0065,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},...
    mfilename,'LR'));

parse(ips,Pi,Hi,Ti,varargin{:});
Ho=ips.Results.Ho;
LR=ips.Results.LR;
clear ips varargin

%% Adjust air pressure
R=287.0; % Ideal gass constant J/kg*K
g=9.81; % Gravitational acceleration m/s2

Pi=readCls(Pi);
Hi=readCls(Hi);
Ti=readCls(Ti);
Ho=readCls(Ho);
LR=readCls(LR);
LR=imresize(LR,size(Pi),'bilinear');

dH=Ho-Hi;
clear Ho Hi
Po=Pi./exp(g*dH./(R*(Ti+LR.*dH/2)));
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end
