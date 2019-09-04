% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 05/10/2019

%% Functionality
% This function calculates two of the three quantities - dew point temperature,
%  specific and relative humidity - given the air temperature and pressure and
%  one of the three.

%% Input
% Ta_fd: details of file or workspace variable for air temperature (K);
% Pa_fd: details of file or workspace variable for air pressure (Pa);
%  InS : characters specifying the type of the third input (Possible types are
%        'Dew point', 'Specific', or 'Relative', in K, g/g, or %);
% In_fd: details of file or workspace variable for the third input.

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
% Require read2Dvar.m and Magnus_F.m.

function [out1,out2]=Cal_Tdw(Ta_fd,Pa_fd,InS,In_fd)
%% Check the inputs
narginchk(4,4);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 4 required inputs\n',mfilename);

addRequired(ips,'Ta_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Ta_fd',1));
addRequired(ips,'Pa_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Pa_fd',2));
expInS={'Specific','Relative','Dew Point'};
msg=cell2mat(cellfun(@(x) [x ', '],expInS,'UniformOutput',false));
msg=sprintf('Expected InS to be one of the following %s\n',msg);
addRequired(ips,'InS',@(x) assert(any(strcmp(x,expInS)),msg));
addRequired(ips,'In_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'In_fd',4));
parse(ips,Ta_fd,Pa_fd,InS,In_fd);
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
Ta=read2Dvar(Ta_fd);
Pa=read2Dvar(Pa_fd);
in3=read2Dvar(In_fd);
clear Ta_fd Pa_fd In_fd
k=isnan(Ta) | isnan(Pa) | isnan(in3);
Ta(k)=NaN;
Pa(k)=NaN;
in3(k)=NaN;
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
es=Magnus_F(Ta);
ms=epsi*es./(Pa-es); % saturated mixing ratio
clear es Aw Ai Bw Bi Cw Ci

switch InS
  case 'Specific' % If specific humidity is known
    e=in3.*Pa./(epsi+(1-epsi)*in3); % from q=epsi*e/(Pa-(1-epsi)*e), in3 is q
    Td=C.*log(e./A)./(B-log(e./A))-abs0; % from Magnus formula
    clear e Pa A B C
    m=in3./(1-in3); % mixing ratio
    RH=m./ms*100; % Relative humidity
    clear m ms in3
    RH(RH>100)=100;

    out1=Td;
    out2=RH;

  case 'Relative' % If relative humidity is known
    m=in3.*ms; % mixing ratio, in3 is RH
    e=m.*Pa./(m+epsi); % from m=epsi*e/(Pa-e)
    Td=C.*log(e./A)./(B-log(e./A))-abs0;
    clear e Pa ms in3 A B C
    q=m./(m+1); % Specific humidity
    clear m

    out1=Td;
    out2=q;

  case 'Dew Point'
    e=Magnus_F(in3); % in3 is Td
    q=epsi*e./(Pa-(1-epsi)*e);
    clear e in3 Pa A B C
    m=q./(1-q); % mixing ratio
    RH=m./ms*100; % Relative humidity
    clear m ms
    RH(RH>100)=100;

    out1=RH;
    out2=q;
end
end
