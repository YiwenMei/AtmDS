% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/2/2019

%% Functionality
% This function is used for the downscaling of incident longwave radiation.

%% Input
% LG_fd : details of file or workspace variable for original incident longwave
%         radiation (W/m2);
% Ta_fd : details of file or workspace variable for original air temperature (K);
% Td_fd : details of file or workspace variable for original dew point temperature (K);
% Tad_fd: details of file or workspace variable for downscaled air temperature (K);
% Tdd_fd: details of file or workspace variable for downscaled dew point temperature (K);
% EmS_cl: characters specifying the methods to calculate atmospheric emissivity
%         (it can be a user-specified emissivity or emissivity calculated based
%         on the Brut, Konz, Satt, Idso, Izio, or Prat method summarized in
%         Fiddes & Grubler (2014). Possible types are 'user', 'brut', 'konz',
%         'satt','idso','izio', or 'prat');

% Emd_fd: if EmS_cl is 'user', use this optional input to specify emissivity
%         as details of file or workspace variable;

%% Output
% LGd: downscaled incident longwave radiation (W/m2);

%% Additional note
% Require read2Dvar.m and Magnus_F.m.

function LGd=LW_DS(LG_fd,Ta_fd,Td_fd,Tad_fd,Tdd_fd,EmS_cl,varargin)
%% Check the inputs
narginchk(6,7);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 6 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'LG_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'LG_fd',1));
addRequired(ips,'Ta_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Ta_fd',2));
addRequired(ips,'Td_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Td_fd',3));
addRequired(ips,'Tad_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Tad_fd',4));
addRequired(ips,'Tdd_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Tdd_fd',5));
addRequired(ips,'EmS_cl',@(x) any(strcmp(x,{'brut','konz','satt','idso','izio','prat'})));

addOptional(ips,'Emd_fd',[],@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Emd_fd',7));
parse(ips,LG_fd,Ta_fd,Td_fd,Tad_fd,Tdd_fd,EmS_cl,varargin{:});
Emd_fd=ips.Results.Emd_fd;
clear ips varargin

%% Read the inputs
FunName=dbstack;
LG=read2Dvar(LG_fd,FunName);
Ta=read2Dvar(Ta_fd,FunName);
Tad=read2Dvar(Tad_fd,FunName);
clear Ta_fd LG_fd Tad_fd

%% Emissivity
switch EmS_cl
  case 'user'
    emd=read2Dvar(Emd_fd,FunName);
    em=imresize(emd,size(LG),'bilinear');
    clear Emd_fd

  otherwise % Calculate emissivity following Fiddes & Grubler (2014)
    sigma=5.670374419e-8; % Stefan-Boltzmann constant (W/m2*K4)
    e=Magnus_F(Td_fd);
    ed=Magnus_F(Tdd_fd);

    switch EmS_cl % Forms of clear-sky emissivity listed in Gubler et al. (2012)
      case 'brut'
        em_cl=1.24*(e./Ta).^(1./7);
        emd_cl=1.24*(ed./Tad).^(1./7);

      case 'konz' % Method addopted in Fiddes & Grubler (2014)
        em_cl=.23+.484*(e./Ta).^(1./8);
        emd_cl=.23+.484*(ed./Tad).^(1./8);

      case 'satt' % Method addopted in Cosgrove et al. (2003)
        em_cl=1.08*(1-exp(-e.^(Ta/2016)));
        emd_cl=1.08*(1-exp(-ed.^(Tad/2016)));

      case 'idso'
        em_cl=.7+5.95e-5*e.*exp(1500./Ta);
        emd_cl=.7+5.95e-5*e.*exp(1500./Tad);

      case 'izio' % Method addopted in Gupta & Tarboton et al. (2016)
        em_cl=1-.43*exp(-11.5*e./Ta);
        emd_cl=1-.43*exp(-11.5*ed./Tad);

      case 'prat'
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
