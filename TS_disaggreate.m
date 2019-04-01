
function tmfl=TS_disaggreate(fTS,fSF,fdn,ndv,cor,wkpth)
% Coarse resolution precipitation
fd_c=read2Dvar(fSF,fdn,ndv);
[~,nm,~]=fileparts(fSF);
% fd_c(fd_c==ndv)=NaN;

% Weighting factor
x=cor(1,2)+cor(5,2)/2:cor(5,2):cor(2,2)-cor(5,2)/2; % Coarse resolution grids
y=cor(3,2)-cor(6,2)/2:-cor(6,2):cor(4,2)+cor(6,2)/2;

fvc=zeros(length(x)*length(y),size(fTS,1)); % Cumulative precipitation
parfor i=1:size(fTS,1)
%   fv=imresize(double(imread(fTS(i,:))),size(fd_c),'bilinear');
  fv=double(imread(fTS(i,:)));
  fvc(:,i)=reshape(fv,numel(fv),1);
end
fv=sum(fvc,2);
k=fv==0;

fvc=fvc./repmat(fv,1,size(fTS,1)); % Weighting factor for precipitation, wp
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
% x=reshape(x,numel(x),1);
% y=reshape(y,numel(y),1);

% [id,d]=knnsearch([x y],[X(k1) Y(k1)],'K',4);
% Id=knnsearch([X Y],[x y]);
% id=Id(id);
[id,d]=knnsearch([x y],[X Y],'K',4);
k1=fd_c>0;
k1=reshape(k1,numel(k1),1);
id=id(k1,:);
d=d(k1,:);
clear x y Id k X Y
d=2.^-(d/35356); % Inversed distance; 25000*sqrt(2)
d=d./sum(d,2); % Weighting factor of wp for the K neighbors

% Time disaggreateion using inversed distance weighted wp
k=reshape(~isnan(fd_c),numel(fd_c),1);
tmfl=cell(size(fTS,1),1);
parfor i=1:size(fTS,1)
  tmfn=[wkpth sprintf('%s%02i.mat',nm,i-1)];
  wp=matfile(tmfn);
  wp=wp.wp;
  wp=sum(wp(id).*d,2);

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
