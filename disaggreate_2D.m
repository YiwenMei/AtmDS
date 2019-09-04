
function Prd=disaggreate_2D(p,Xc,Yc,pd,Xdc,Ydc)
%% Associate DS grids to target grids
p=read2Dvar(p);
Xc=read2Dvar(Xc);
Xc=Xc(~isnan(p));
% Xc=reshape(Xc,numel(Xc),1);
Yc=read2Dvar(Yc);
% Yc=reshape(Yc,numel(Yc),1);
Yc=Yc(~isnan(p));

pd=read2Dvar(pd);
Xdc=read2Dvar(Xdc);
Xdc=Xdc(~isnan(pd));
Ydc=read2Dvar(Ydc);
Ydc=Ydc(~isnan(pd));

id=knnsearch([Xc Yc],[Xdc Ydc]);
clear Xc Yc

%% Find the weighting factor
pdu=accumarray(id,pd(~isnan(pd)),[],[],NaN); % Sum DS grids associated to a same target grid
n=accumarray(id,pd(~isnan(pd))>0); % Count DS grids with non-zero rainfall associated to a same target grid
N=accumarray(id,1); % Count DS grids associated to a same target grid
pdu=pdu./n;
pdu(n==0)=0;

pdu=repelem(pdu,N); % Repeat the elements
[~,I]=sort(id);
pdu(I)=pdu; % Sort the elements back to the order of id

W=pd(~isnan(pd))./pdu;
W(pdu==0)=0;
clear pdu n N I

%% Find the case of miss (CHIRPS>0 but DS=0)
% p=read2Dvar(p);
pf=p(~isnan(p));
pi=pf(id);
% pi=p(id);
id=knnsearch([Xdc(W>0) Ydc(W>0)],[Xdc(pi>0 & W==0) Ydc(pi>0 & W==0)]); % Closest DS grid with p>0

w=W(W>0);
W(pi>0 & W==0)=w(id);
w=pd;
w(~isnan(pd))=W; % Map the weighting factor to the target grids
clear id W

%% Apply the weighting factor to CHIRPS
pk=pd;
pk(~isnan(pd))=pi;
Prd=w.*pk;
end
