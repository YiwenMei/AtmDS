% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/16/2019

%% Functionality
% This function calculates two of the three quantities - dew point temperature,
%  specific and relative humidity - given the air temperature and pressure and
%  one of the three.

%% Input
% Ta : spatial map class (V2DCls.m) object or workspace variable for air temperature
%       (K);
% Pa : V2DCls.m object or workspace variable for air pressure (Pa);
% InS: characters specifying the type of the third input (Possible types are
%       'Dew point', 'Specific', or 'Relative', in K, g/g, or %);
% InV: V2DCls.m object or workspace variable for the third input.

%% Output
% If specific humidity is supplied
% out1: dew point temperature (K);
% out2: relative humidity (%);

% If relative humidity is supplied
% out1: dew point temperature (K);
% out2: specific humidity (g/g);

% If dew point temperature is supplied
% out1: relative humidity (%);
% out2: specific humidity (g/g).

%% Additional note
% Require V2DCls.m.

function [out1,out2]=Hum_Cal(Ta,Pa,InS,InV)
%% Check the inputs
narginchk(4,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Ta',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Ta'));
addRequired(ips,'Pa',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Pa'));
expInS={'Specific','Relative','Dew Point'};
msg=cell2mat(cellfun(@(x) [x ', '],expInS,'UniformOutput',false));
msg=sprintf('Expected InS to be one of the following %s\n',msg);
addRequired(ips,'InS',@(x) assert(any(strcmp(x,expInS)),msg));
addRequired(ips,'InV',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'InV'));

parse(ips,Ta,Pa,InS,InV);
clear ips msg

%% Constant
epsi=.62198; % Ratio of molecular weight of water and dry air
abs0=-273.15;
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

%% Read the inputs
Ta=readCls(Ta);
Pa=readCls(Pa);
InV=readCls(InV);
k=isnan(Ta) | isnan(Pa) | isnan(InV);
Ta(k)=NaN;
Pa(k)=NaN;
InV(k)=NaN;
clear k

%% Calculate Tdw
A=Am*ones(size(Ta));
A(Ta>5)=Aw;
A(Ta<-5)=Ai;
B=Bm*ones(size(Ta));
B(Ta>5)=Bw;
B(Ta<-5)=Bi;
C=Cm*ones(size(Ta));
C(Ta>5)=Cw;
C(Ta<-5)=Ci;
es=A.*exp(B.*(Ta+abs0)./(Ta+abs0+C)); % Magnus formula
ms=epsi*es./(Pa-es); % saturated mixing ratio
clear es Aw Ai Bw Bi Cw Ci

switch InS
  case 'Specific' % If specific humidity is known
    e=InV.*Pa./(epsi+(1-epsi)*InV); % from q=epsi*e/(Pa-(1-epsi)*e), InV is q
    Td=C.*log(e./A)./(B-log(e./A))-abs0; % from Magnus formula
    clear e Pa A B C
    m=InV./(1-InV); % mixing ratio
    RH=m./ms*100; % Relative humidity
    clear m ms InV
    RH(RH>100)=100;

    out1=Td;
    out2=RH;

  case 'Relative' % If relative humidity is known
    m=InV.*ms; % mixing ratio, InV is RH
    e=m.*Pa./(m+epsi); % from m=epsi*e/(Pa-e)
    Td=C.*log(e./A)./(B-log(e./A))-abs0;
    clear e Pa ms InV A B C
    q=m./(m+1); % Specific humidity
    clear m

    out1=Td;
    out2=q;

  case 'Dew Point'
    e=A.*exp(B.*(InV+abs0)./(InV+abs0+C)); % InV is Td
    q=epsi*e./(Pa-(1-epsi)*e);
    clear e InV Pa A B C
    m=q./(1-q); % mixing ratio
    RH=m./ms*100; % Relative humidity
    clear m ms
    RH(RH>100)=100;

    out1=RH;
    out2=q;
end
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end
