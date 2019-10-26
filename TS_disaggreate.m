% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 5/25/2019

%% Functionality
% This function disaggreate a spatially high but temporally coarse field 1 by
%  a temporally high but spatially coarse field 2 (e.g. disaggreate 1km/daily
%  precipitation to 1km/hourly by another 50km/hourly precipitation). Note that
%  field 1 needs to be covered entirely by field 2.

%% Input
%  fTS : a space-time variable class (V2DTCls.m) object for the high temporal
%        resolution fields;
%  fSF : a matlab variable for the high spatial resolution field;
%  cSF : boundary of the high spatial resolution fields (cSF is 3-by-2 as follow
%        [xl yt;xr yb;Rx Ry], where x/y is the horizontal/vertical dimension,
%        l/r/t/b is the left/right/top/bottom boundary, and R is the resolution);
% wkpth: working directory of the code.

% pflg: parallel flag (false - default, squential; true - parallel).

%% Output
% Tfl: a cell array stores the file names for the disaggreated field located
%      in wkpth.

%% Adittional note
% Require V2DTCls.m.

function Tfl=TS_disaggreate(fTS,fSF,cSF,wkpth,varargin)
%% Check the inputs
narginchk(4,5);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'fTS',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'fTS'));
addRequired(ips,'fSF',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'fSF'));
addRequired(ips,'cSF',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'cSF'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,fTS,fSF,cSF,wkpth,varargin{:});
pflg=ips.Results.pflg;
clear ips varargin

%% Hourly weighting factor of precipitation
fvc=[];
for i=1:length(fTS.Fnm)
  fv=fTS.readCls(i);
  fvc=nansum(cat(3,fvc,fv),3);
end
clear fv

Tfl=cell(length(fTS.Fnm),1);
switch pflg
  case true
    parfor i=1:length(fTS.Fnm)
      Tfl{i}=TS_dis_sub1(fTS,i,fvc,wkpth);
    end

  case false
    for i=1:length(fTS.Fnm)
      Tfl{i}=TS_dis_sub1(fTS,i,fvc,wkpth);
    end
end

%% Time disaggregation
[X,Y]=meshgrid(cSF(1,1)+cSF(3,1)/2:cSF(3,1):cSF(2,1)-cSF(3,1)/2,...
  cSF(1,2)-cSF(3,2)/2:-cSF(3,2):cSF(2,2)+cSF(3,2)/2); % High resolution grids
X=X(fSF>0);
Y=Y(fSF>0);

% Coarse resolution precipitation
[x,y]=fTS.GridCls; % Coarse resolution grids
x=x(fvc>0);
y=y(fvc>0);

[id,d]=knnsearch([x y],[X Y],'K',4);
clear x y X Y

d=2.^-(d/hypot(fTS.Gtg(3,1)/2,fTS.Gtg(3,2)/2)); % Inversed distance;
d=d./sum(d,2); % Weighting factor of wp for the K neighbors

% Time disaggreateion using inversed distance weighted wp
switch pflg
  case true
    parfor i=1:length(fTS.Fnm)
      TS_dis_sub2(Tfl{i},wkpth,fvc,id,d,fSF);
    end

  case false
    for i=1:length(fTS.Fnm)
      TS_dis_sub2(Tfl{i},fvc,id,d,fSF);
    end
end
end

function [tmfn,wp]=TS_dis_sub1(fTS,n,fvc,wkpth)
wp=fTS.readCls(n);
wp=wp./fvc; % Hourly weighting factor, wp
wp(fvc==0)=0;

[~,nm,~]=fileparts(fTS.Fnm{n});
tmfn=fullfile(wkpth,sprintf('w%s.mat',nm));
save(tmfn,'wp');
end

function TS_dis_sub2(tfn,fvc,id,d,fSF)
wp=matfile(tfn);
wp=wp.wp;
wp=wp(fvc>0);
wp=wp(id);
wp=sum(wp.*d,2);

pr=fSF;
pr(fSF>0)=wp;
pr=pr.*fSF;

save(tfn,'pr');
end
