% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 03/29/2019

%% Functionality
% This function is used for the downscaling of incident longwave radiation. It
% supports both user-specific emissivity variables or emissivity calculated using
% Cosgrove et al. (2003).

%% Input
% LGfn: full name of file or array for coarse resolution incident longwave radiation
%       at land surface (W/m2);
% Tafn: full name of file or array for coarse resolution air temperature (K);
% Tad : downscaled air temperature (K);
% ndv : no-data value for the inputs dataset (use only one ndv for all inputs);
% pr1 : coarse resolution emissivity or full name of file or array for coarse
%       resolution air pressure (Pa);
% pr2 : downscaled emissivity or downscaled air pressure (Pa);
% pr3 : full name of file or array for coarse resolution specific humidity (g/g);
% pr4 : downscaled specific humidity (g/g);

%% Output
% LGd : downscaled incident longwave radiation at land surface (W/m2);

function LGd=LW_DS(LGfn,Tafn,Tad,ndv,pr1,pr2,pr3,pr4)
%% Check the inputs
epsi=.62198; % Ratio of molecular weight of water and dry air
switch nargin
    case {1:5}; error('Not enough arguments');
    case 6
        EmT='user';
        em=pr1;
        emd=pr2;
        clear pr1 pr2
    case 7; error('Not enough arguments for C2003 method');
    case 8
        EmT='C2003';
        Pafn=pr1;
        Pad=pr2;
        qfn=pr3;
        qd=pr4;
        clear pr1 pr2 pr3 pr4
    otherwise; error('Too many number of arguments');
end

%% Read the inputs
if ischar(Tafn)
  Ta=double(imread(Tafn));
else
  Ta=Tafn;
end
Ta(Ta==ndv)=NaN;
if ischar(LGfn)
  LG=double(imread(LGfn));
else
  LG=LGfn;
end
LG(LG==ndv)=NaN;
clear Tafn LGfn Tafn Pafn qfn

%% C2003 emissivity
if strcmp(EmT,'C2003')
  if ischar(Pafn)
    Pa=double(imread(Pafn));
  else
    Pa=Pafn;
  end
  Pa(Pa==ndv)=NaN;
  if ischar(qfn)
    q=double(imread(qfn));
  else
    q=qfn;
  end
  q(q==ndv)=NaN;

  ed=qd.*Pad./(epsi+(1-epsi)*pr4);
  emd=1.08*(1-exp(-ed.^(Tad/2016))); % Cosgrove et al. (2003) emissivity
  e=q.*Pa./(epsi+(1-epsi)*q);
  em=1.08*(1-exp(-e.^(Ta/2016))); % Cosgrove et al. (2003) emissivity
  clear ed qd Pad e q Pa
end

%% Downscaling of longwave radiation
kem=emd./imresize(em,size(emd),'bilinear');
te=log(kem);
te1=te(~isnan(te));
tem=mean(te1);
tu=quantile(te(te>tem),.995);
td=quantile(te(te<tem),.005);
kem(te>tu)=exp(tu);
kem(te<td)=exp(td);
kTa=Tad./imresize(Ta,size(Tad),'bilinear'); % LG=ems*sigma*Ta^4
LGd=kem.*kTa.^4.*imresize(LG,size(Tad),'bilinear');
end
