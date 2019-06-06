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
%  fTS : full name of files for the high temporal resolution precipitation fields;
%  fSF : full name of file for the high spatial resolution precipitation field;
%  fdn : field name of fSF if fSF is .nc, .hdf or .mat file;
%  ndv : no-data value for the inputs dataset (use only one ndv for all inputs);
%  cor : coordinates of the two precipitation fields (cor is 6-by-2 where 6 represents
%        the the left, right, top, bottom, x resolution, and y resolution and
%        2 represents the coordinates of fSF and fTS);
%   K  : the K neighbors contributing to the target grid cell;
% wkpth: working directory of the code.

%% Output
% tmfl: full name of files for the outputted high spatiotemporal precipitation
%       fields in .mat.

function tmfl=TS_disaggreate(fTS,fSF,fdn,ndv,cor,K,wkpth)
% Coarse resolution precipitation
fd_c=read2Dvar(fSF,fdn,ndv);
[~,nm,~]=fileparts(fSF);

x=cor(1,2)+cor(5,2)/2:cor(5,2):cor(2,2)-cor(5,2)/2; % Coarse resolution grids
y=cor(3,2)-cor(6,2)/2:-cor(6,2):cor(4,2)+cor(6,2)/2;

% Hourly weighting factor of precipitation
fvc=zeros(length(x)*length(y),size(fTS,1));
parfor i=1:size(fTS,1)
  fv=double(imread(fTS(i,:)));
  fvc(:,i)=reshape(fv,numel(fv),1);
end
fv=sum(fvc,2); % Cumulative precipitation
k=fv==0;

fvc=fvc./repmat(fv,1,size(fTS,1)); % Hourly weighting factor, wp
parfor i=1:size(fTS,1)
  tmfn=[wkpth sprintf('%s%02i.mat',nm,i-1)];
  parsave(tmfn,fvc(:,i),'wp');
end
clear fvc fv

%% Time disaggregation
[X,Y]=meshgrid(cor(1,1)+cor(5,1)/2:cor(5,1):cor(2,1)-cor(5,1)/2,...
    cor(3,1)-cor(6,1)/2:-cor(6,1):cor(4,1)+cor(6,1)/2);
X=reshape(X,numel(X),1); % High resolution grids
Y=reshape(Y,numel(Y),1);

[x,y]=meshgrid(x,y); % Coarse resolution grids
x=x(~k);
y=y(~k);

[id,d]=knnsearch([x y],[X Y],'K',K);
k1=fd_c>0;
k1=reshape(k1,numel(k1),1);
id=id(k1,:);
d=d(k1,:);
clear x y Id k X Y
d=2.^-(d/hypot(cor(5,2)/2,cor(6,2)/2)); % Inversed distance;
d=d./sum(d,2); % Weighting factor of wp for the K neighbors

% Time disaggreateion using inversed distance weighted wp
k=reshape(~isnan(fd_c),numel(fd_c),1);
tmfl=cell(size(fTS,1),1);
parfor i=1:size(fTS,1)
  tmfn=[wkpth sprintf('%s%02i.mat',nm,i-1)];
  wp=matfile(tmfn);
  wp=wp.wp;
  wp=wp(id);
  K=sum(~isnan(wp),2);
  wd=d;
  wd(isnan(wp))=0;
  wd=wd./sum(wd,2); % Updated weighting factor of wp
  wp=nanmean(wp.*wd,2).*K;
  wp(K==0)=1/size(fTS,1);

  Wp=nan(size(k1));
  Wp(k1)=wp;
  Wp(~k1 & k)=0;
  Wp=reshape(Wp,size(fd_c));
  fd_h=fd_c.*Wp;
  fd_h(isnan(fd_h))=ndv;

% Save the files
  parsave(tmfn,fd_h,'pr');
  tmfl{i}=tmfn;
end
tmfl=cell2mat(tmfl);
end
