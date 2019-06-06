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

function [LGd,emd,sts]=LW_DS(LGfn,Tafn,Tad,ndv,pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10,pr11)
%% Check the inputs
epsi=.62198; % Ratio of molecular weight of water and dry air
sigma=5.670374419e-8; % Stefan-Boltzmann constant (W/m2*K4)
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
    case {9:14}; error('Not enough arguments for A2007 method');
    case 15
        EmT='A2007';
        Pafn=pr1;
        Pad=pr2;
        SGfn=pr3;
        STfn=pr4;
        Tdfn=pr5;
        Tdd=pr6;
        El=pr7;
        Eld=pr8;
        Azd=pr9;
        Asp=pr10;
        Slp=pr11;
        clear pr1 pr2 pr3 pr4 pr5 pr6 pr7 pr8 pr9 pr10 pr11;
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
clear Tafn LGfn
sts=[];

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
  clear Pafn qfn

  ed=qd.*Pad./(epsi+(1-epsi)*qd);
  emd=1.08*(1-exp(-ed.^(Tad/2016))); % Cosgrove et al. (2003) emissivity
  e=q.*Pa./(epsi+(1-epsi)*q);
  em=1.08*(1-exp(-e.^(Ta/2016))); % Cosgrove et al. (2003) emissivity
  clear ed qd Pad e q Pa

  sts=[length(find(emd>1))/length(find(~isnan(emd))) length(find(em>1))/length(find(~isnan(em)))];
end

%% A2007 emissivity
if strcmp(EmT,'A2007')
  if ischar(Pafn)
    Pa=double(imread(Pafn));
  else
    Pa=Pafn;
  end
  Pa(Pa==ndv)=NaN;
  if ischar(SGfn)
    SG=double(imread(SGfn));
  else
    SG=SGfn;
  end
  SG(SG==ndv)=NaN;
  if ischar(STfn)
    ST=double(imread(STfn));
  else
    ST=STfn;
  end
  ST(ST==ndv)=NaN;
  if ischar(Tdfn)
    Td=double(imread(Tdfn));
  else
    Td=Tdfn;
  end
  Td(Td==ndv)=NaN;
  clear Pafn SGfn STfn Tdfn

  Asp(Asp==ndv)=NaN;
  Slp(Slp==ndv)=NaN;
  El(El<0)=0;
  El(El>90)=90;
  ze=90-El;
  ze(ze>85)=85;
  Eld(Eld==ndv)=NaN;
  Eld(Eld<0)=0;
  Eld(Eld>90)=90;
  zed=90-Eld;
  zed(zed>85)=85;
  Azd(Azd==ndv)=NaN;
  clear El Eld

% Coarse resolution broad-band atmospheric transmissivity
  kt=SG./ST;
  cosi=cosd(ze);
  te=-.00146*(Pa/1000)./cosi./kt;
  e=Magnus_F(Td)/1000;
  W=.14*e.*Pa/1000+2.1;
  te1=-.075*(W./cosi).^.4;
  tau_sw=.35+.627*exp(te+te1); % Broad-band atmospheric transmissivity
  tau_sw=reshape(tau_sw,numel(tau_sw),1);

% Model fitting
  em=LG./Ta.^4/sigma;
  y=log(reshape(em,numel(em),1));
  X=[ones(size(y)) log(-log(tau_sw))];
  [b,~,res,~,sts]=regress(y,X); % Allen et al. (2007) Eq.(25)
  sts=[b' sqrt(sum(res.^2)/(length(res)-length(b))) sts];

% High resolution broad-band atmospheric transmissivity
  cosid=cosd(Slp).*cosd(zed)+sind(Slp).*sind(zed).*cosd(Azd-Asp);
  cosid(Slp==0)=cosd(zed(Slp==0)); % Slp=0 means Asp=NaN except the ocean
  cosid(cosid<cosd(85))=cosd(85);
  te=-.00146*(Pad/1000)./cosid./imresize(kt,size(Pad),'bilinear');
  e=Magnus_F(Tdd)/1000;
  W=.14*e.*Pad/1000+2.1;
  te1=-.075*(W./cosid).^.4;
  tau_sw=.35+.627*exp(te+te1); % Broad-band atmospheric transmissivity
  clear cosid Slp zen Azd Asp e W te te1

% Estimate high resolution emissivity
  tau_sw=reshape(tau_sw,numel(tau_sw),1);
  y=[ones(size(tau_sw)) log(-log(tau_sw))]*b;
  emd=exp(reshape(y,size(Tad)));
end

%% Downscaling of longwave radiation
kem=emd./imresize(em,size(emd),'bilinear');
te=log(kem);
te1=te(~isnan(te));
thr=1:.1:6;
pfg=0; % Turn to 1 for a check
[~,tu,td]=RemOut_2D(te1,thr,pfg); 
kem(te>tu)=exp(tu);
kem(te<td)=exp(td);
kTa=Tad./imresize(Ta,size(Tad),'bilinear'); % LG=ems*sigma*Ta^4
LGd=kem.*kTa.^4.*imresize(LG,size(Tad),'bilinear');
end
