
function [Em,Emd,sts]=Cal_emR(LG,Ta,SG,ST,SGd,Sdd)
a=.952; % Coefficient from Ruiz-Arias et al. (2012)
b=1.041;
c=2.3;
d=4.702;
l=.952-1.041*exp(-exp(2.3-4.702*1));
m=.952-1.041*exp(-exp(2.3-4.702*0));
sigma=5.67036713e-8; % Stefan-Boltzmann constant, W/(m^2*K^4) 

% Estimate the model parameters
LGa=mean(LG,3);
Taa=mean(Ta,3);
Em=LGa./(sigma*(Taa).^4); % MERRA2 emissivity derived by Stefan-Boltzmann law
y=reshape(log(Em),numel(Em),1);
kt=SG./ST; % Clear sky index
X=[reshape(log(-log(kt)),numel(kt),1) ones(numel(kt),1)];
[h,~,res,~,sts]=regress(y,X); % Allen et al. (2007) model
clear SG ST Taa LGa

% High resolution emissivity
res=reshape(res,size(Em)); % Residual
kd=Sdd./SGd; % Weighting factor of diffuse shortwave
kd=(m-l)*kd+l;
kt=c/d-log(log(1./(a/b-kd/b)))/d; % Inversed Ruiz-Arias et al. (2012) regression model
X=[reshape(log(-log(kt)),numel(kt),1) ones(numel(kt),1)];
y=X*h;
y=reshape(y,size(SGd));
Emd=exp(y+imresize(res,size(SGd),'bilinear')); % High resolution emissivity
end
