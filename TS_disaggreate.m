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

function tmfl=TS_disaggreate(fTS,cor,fSF,cor1,K,wkpth)
% Hourly weighting factor of precipitation
fvc=[];
for i=1:size(fTS{1},1)
  fv=fTS;
  fv{1}=fTS{1}(i,:);
  fv=read2Dvar(fv);
  fvc=nansum(cat(3,fvc,fv),3);
end

parfor i=1:size(fTS{1},1)
  fv=fTS;
  fv{1}=fTS{1}(i,:);
  fv=read2Dvar(fv);
  wp=fv./fvc; % Hourly weighting factor, wp
  wp(fvc==0)=0;

  [~,nm,~]=fileparts(fTS{1}(i,:));
  tmfn=[wkpth sprintf('w%s.mat',nm)];
  parsave(tmfn,wp,'wp');
end
clear fv wp

%% Time disaggregation
[X,Y]=meshgrid(cor1(1,1)+cor1(3,1)/2:cor1(3,1):cor1(2,1)-cor1(3,1)/2,...
    cor1(1,2)-cor1(3,2)/2:-cor1(3,2):cor1(2,2)+cor1(3,2)/2);
fSF=read2Dvar(fSF);
X=X(fSF>0);
Y=Y(fSF>0);
% X=reshape(X,numel(X),1); % High resolution grids
% Y=reshape(Y,numel(Y),1);

% Coarse resolution precipitation
x=cor(1,1)+cor(3,1)/2:cor(3,1):cor(2,1)-cor(3,1)/2; % Coarse resolution grids
y=cor(1,2)-cor(3,2)/2:-cor(3,2):cor(2,2)+cor(3,2)/2;
[x,y]=meshgrid(x,y); % Coarse resolution grids
x=x(fvc>0);
y=y(fvc>0);

[id,d]=knnsearch([x y],[X Y],'K',K);
clear x y X Y

d=2.^-(d/hypot(cor(3,1)/2,cor(3,2)/2)); % Inversed distance;
d=d./sum(d,2); % Weighting factor of wp for the K neighbors

% Time disaggreateion using inversed distance weighted wp
tmfl=cell(size(fTS{1},1),1);
parfor i=1:size(fTS{1},1)
  [~,nm,~]=fileparts(fTS{1}(i,:));
  tmfn=[wkpth sprintf('w%s.mat',nm)];
  wp=matfile(tmfn);
  wp=wp.wp;
  wp=wp(fvc>0);
  wp=wp(id);
  wp=sum(wp.*d,2);

  fd_h=fSF;
  fd_h(fSF>0)=wp;
  fd_h=fd_h.*fSF;
  fd_h(isnan(fd_h))=fTS{2};

% Save the files
  parsave(tmfn,fd_h,'pr');
  tmfl{i}=tmfn;
end
tmfl=cell2mat(tmfl);
end
