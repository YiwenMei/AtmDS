% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 5/25/2019

%% Functionality
% This function disaggreate spatially high resolution cumulative precipitation
% by temporally high resolution precipitation fields (e.g. disaggreate 1km/daily
% precipitation to 1km/hourly by another 50km/hourly precipitation fields). Note
% that the temporally high resolution fields need to cover a larger extent than
% the spatial one.

%% Input
%  fTS : file name list with details of the high temporal resolution fields;
%  cTS : boundary of the high temporal resolution fields (cTS is 3-by-2 as follow
%        [xl yt;xr yb;Rx Ry], where x/y is the horizontal/vertical dimension,
%        l/r/t/b is the left/right/top/bottom boundary, and R is the resolution);
%  fSF : details of file or array for the high spatial resolution field;
%  cSF : same as cTS but for SF;
%   K  : the K neighbors from TS grid contributing to the SF grid;
% wkpth: working directory of the code.

%% Output
% tmfl: file name list for the outputted high spatiotemporal fields in .mat.

%% Adittional note
% Require read2Dvar.m and parsave.m.

function tmfl=TS_disaggreate(fTS,cTS,fSF,cSF,K,wkpth,pflg)
%% Hourly weighting factor of precipitation
fvc=[];
for i=1:size(fTS{1},1)
  fv=fTS;
  fv{1}=fTS{1}(i,:);
  fv=read2Dvar(fv);
  fvc=nansum(cat(3,fvc,fv),3);
end

switch pflg
  case true
    parfor i=1:size(fTS{1},1)
      fv=fTS;
      fv{1}=fTS{1}(i,:);
      [~]=TS_dis_sub1(fv,fvc,wkpth);
    end

  case false
    for i=1:size(fTS{1},1)
      fv=fTS;
      fv{1}=fTS{1}(i,:);
      [~]=TS_dis_sub1(fv,fvc,wkpth);
    end
end
clear fv

%% Time disaggregation
[X,Y]=meshgrid(cSF(1,1)+cSF(3,1)/2:cSF(3,1):cSF(2,1)-cSF(3,1)/2,...
    cSF(1,2)-cSF(3,2)/2:-cSF(3,2):cSF(2,2)+cSF(3,2)/2); % High resolution grids
fSF=read2Dvar(fSF);
X=X(fSF>0);
Y=Y(fSF>0);

% Coarse resolution precipitation
[x,y]=meshgrid(cTS(1,1)+cTS(3,1)/2:cTS(3,1):cTS(2,1)-cTS(3,1)/2,...
    cTS(1,2)-cTS(3,2)/2:-cTS(3,2):cTS(2,2)+cTS(3,2)/2); % Coarse resolution grids
x=x(fvc>0);
y=y(fvc>0);

[id,d]=knnsearch([x y],[X Y],'K',K);
clear x y X Y

d=2.^-(d/hypot(cTS(3,1)/2,cTS(3,2)/2)); % Inversed distance;
d=d./sum(d,2); % Weighting factor of wp for the K neighbors

% Time disaggreateion using inversed distance weighted wp
tmfl=cell(size(fTS{1},1),1);
switch pflg
  case true
    parfor i=1:size(fTS{1},1)
      tmfn=TS_dis_sub2(fTS{1}(i,:),wkpth,fvc,id,d,fSF,fTS{2});
      tmfl{i}=tmfn;
    end

  case false
    for i=1:size(fTS{1},1)
      tmfn=TS_dis_sub2(fTS{1}(i,:),wkpth,fvc,id,d,fSF,fTS{2});
      tmfl{i}=tmfn;
    end
end
tmfl=cell2mat(tmfl);
end

function wp=TS_dis_sub1(TS,fvc,wkpth)
wp=read2Dvar(TS);
wp=wp./fvc; % Hourly weighting factor, wp
wp(fvc==0)=0;

[~,nm,~]=fileparts(TS{1});
tmfn=fullfile(wkpth,sprintf('w%s.mat',nm));
save(tmfn,'wp');
end

function tmfn=TS_dis_sub2(fTS,wkpth,fvc,id,d,fSF,ndv)
[~,nm,~]=fileparts(fTS);
tmfn=fullfile(wkpth,sprintf('w%s.mat',nm));
wp=matfile(tmfn);
wp=wp.wp;
wp=wp(fvc>0);
wp=wp(id);
wp=sum(wp.*d,2);

pr=fSF;
pr(fSF>0)=wp;
pr=pr.*fSF;
pr(isnan(pr))=ndv;

% Save the files
save(tmfn,'pr');
end
