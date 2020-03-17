function LGd=LW_DS(LG,Ta,Td,Tad,Tdd,EmcT,varargin)
%% Check the inputs
narginchk(6,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'LG',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'LG'));
addRequired(ips,'Ta',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Ta'));
addRequired(ips,'Td',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Td'));
addRequired(ips,'Tad',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Tad'));
addRequired(ips,'Tdd',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Tdd'));
addRequired(ips,'EcT',@(x) any(strcmp(x,{'User','Brut','Konz','Satt','Idso','Izio','Prat'})));

addOptional(ips,'Emd',[],@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Emd'));
parse(ips,LG,Ta,Td,Tad,Tdd,EmcT,varargin{:});
Emd=ips.Results.Emd;
clear ips varargin

%% Read the inputs
LG=readCls(LG);
Ta=readCls(Ta);
Tad=readCls(Tad);
Tdd=readCls(Tdd);

%% Emissivity
switch EmcT
  case 'User'
    emd=readCls(Emd);
    em=imresize(emd,size(LG),'bilinear');

  otherwise % Calculate emissivity following Fiddes & Grubler (2014)
    sigma=5.670374419e-8; % Stefan-Boltzmann constant (W/m2*K4)
    e=Magnus_F(Td,Ta);
    ed=Magnus_F(Tdd,Tad);

    switch EmcT % Forms of clear-sky emissivity listed in Gubler et al. (2012)
      case 'Brut'
        em_cl=1.24*(e./Ta).^(1./7);
        emd_cl=1.24*(ed./Tad).^(1./7);

      case 'Konz' % Method addopted in Fiddes & Grubler (2014)
        em_cl=.23+.484*(e./Ta).^(1./8);
        emd_cl=.23+.484*(ed./Tad).^(1./8);

      case 'Satt' % Method addopted in Cosgrove et al. (2003)
        em_cl=1.08*(1-exp(-e.^(Ta/2016)));
        emd_cl=1.08*(1-exp(-ed.^(Tad/2016)));

      case 'Idso'
        em_cl=.7+5.95e-5*e.*exp(1500./Ta);
        emd_cl=.7+5.95e-5*e.*exp(1500./Tad);

      case 'Izio' % Method addopted in Gupta & Tarboton et al. (2016)
        em_cl=1-.43*exp(-11.5*e./Ta);
        emd_cl=1-.43*exp(-11.5*ed./Tad);

      case 'Prat'
        em_cl=1-(1+46.5*e./Ta)*exp(-(1.2+3*46.5*e./Ta).^.5);
        emd_cl=1-(1+46.5*ed./Tad)*exp(-(1.2+3*46.5*ed./Tad).^.5);
    end

% Additive All-sky emissivity
    em=LG./Ta.^4/sigma; % All-sky emissivity
    em_c=em-em_cl; % Cloud emissivity
    emd=emd_cl+imresize(em_c,size(Tad),'bilinear');
    if ~isempty(find(emd<0, 1))
      emi=imresize(em,size(Tad),'bilinear');
      emd(emd<0)=emi(emd<0);
    end
    clear e ed em_c emi
end

%% Downscaling of longwave radiation
kem=emd./imresize(em,size(emd),'bilinear');
kTa=Tad./imresize(Ta,size(Tad),'bilinear'); % LG=ems*sigma*Ta^4
LGd=kem.*kTa.^4.*imresize(LG,size(Tad),'bilinear');
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

function E=Magnus_F(T,Ta)
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

A=Am*ones(size(Ta));
A(Ta>5)=Aw;
A(Ta<-5)=Ai;
B=Bm*ones(size(Ta));
B(Ta>5)=Bw;
B(Ta<-5)=Bi;
C=Cm*ones(size(Ta));
C(Ta>5)=Cw;
C(Ta<-5)=Ci;

E=A.*exp(B.*(T+abs0)./(T+abs0+C)); % Magnus formula
end
