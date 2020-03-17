function VI=Tinterp2D(imb,ima,dnb,dna,dn)
%% Check inputs
narginchk(5,5);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'imb',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'imb'));
addRequired(ips,'ima',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'ima'));
addRequired(ips,'dnb',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dnb'));
addRequired(ips,'dna',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dna'));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dn'));

parse(ips,imb,ima,dnb,dna,dn);
clear ips

%% Mean VI
imb=readCls(imb);
ima=readCls(ima);
a=size(ima);

imb(isnan(imb))=ima(isnan(imb)); % Fill NaN by each other's records
ima(isnan(ima))=imb(isnan(ima));
imb=reshape(imb,[numel(imb) 1]);
ima=reshape(ima,[numel(ima) 1]);

w=1-abs(dn-[dnb;dna])/diff([dnb;dna]);
VI=reshape([imb ima]*w,a);
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
