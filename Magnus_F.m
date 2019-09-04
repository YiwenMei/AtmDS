function E=Magnus_F(T_fd,varargin)
%% Check inputs
narginchk(1,2);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'T_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'T_fd',1));
addOptional(ips,'mflg','ABC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'mflg',2));
parse(ips,T_fd,varargin{:});
mflg=ips.Results.mflg;
clear ips

%% Inputs
FunName=dbstack;
T=read2Dvar(T_fd,FunName);

%% Calculation
abs0=-273.15; % K
switch mflg

  case 'ABC'
% Coefficient for Magnus formula adpoted from Buck (1981)
    Aw=611.21; % If Ta>5, use coeffcient of the ew2 curve in Buck (1981)
    Bw=17.368;
    Cw=238.88;
    Am=611.21; % If -5<=Ta<=5, use coeffcient of the ew1 curve in Buck (1981)
    Bm=17.502;
    Cm=240.97;
    Ai=611.15; % If Ta<-5, use coeffcient of the ei2 curve in Buck (1981)
    Bi=22.452;
    Ci=272.55;

% Saturated vapor pressure (Pa)
    E=Aw*exp(Bw*(T+abs0)./(T+abs0+Cw)); % Magnus formula
    E(T<=5-abs0)=Am*exp(Bm*(T(T<=5-abs0)+abs0)./(T(T<=5-abs0)+abs0+Cm));
    E(T<-5-abs0)=Ai*exp(Bi*(T(T<-5-abs0)+abs0)./(T(T<-5-abs0)+abs0+Ci));

  case 'L&R'
    R=461.5; % J/K/kg
    es0=exp(.0511*T+6.686)/100; % Pa
    L=2500800-2360*T+1.6*T.^2-.06*T.^3; % J/kg
    L(T<0)=2834100-290*T(T<0)-4*T(T<0).^2;

% Saturated vapor pressure (Pa)
    E=es0.*exp((L./R)*(1/-abs0-1./T)); % Magnus formula
end
end
