function AtmDS_exe(inifn,pflg)
%% Read the control file
fprintf('Read the control file\n');
ini=IniConfig();
ini.ReadFile(inifn);
Sec=ini.GetSections();
for i=1:length(Sec)
  keys=ini.GetKeys(Sec{i});
  value=ini.GetValues(Sec{i},keys);
  for j=1:length(keys)
    if i==1 && j<=2
      cmdstr=sprintf('%s.%s=value{j}(2:end);',Sec{i}(2:end-1),keys{j});
    else
      cmdstr=sprintf('%s.%s=value{j};',Sec{i}(2:end-1),keys{j});
    end
    eval(cmdstr);
  end
end
clear i j keys value Sec ini cmdstr inifn

%% Run option
fprintf('Selected model options:\n');
if isscalar(Forcing.MHti_Vb) && Forcing.MHti_Vb==Forcing.MHto_Vb % Adjust measurement height or not
  fprintf('1. Input and output height are the same\n');
  flg=false;
else
  fprintf('1. Adjust measurement height to %.1fm\n',Forcing.MHto_Vb);
  flg=true;
end

if strcmp(Parameter.LRTp_Vb,'User')
  fprintf('2. Air temperature lapse rate provided\n');
else
  fprintf('2. Calculate air temperature lapse rate by %s method\n',Parameter.LRTp_Vb);
end

if strcmp(Forcing.IVTp_Vb,'DewPoint') % Dew point, specific or relative humidity
  fprintf('3. Dew point temperature provided\n');
  if strcmp(Parameter.LRTp_Vb,'User')
    fprintf(' - Dew point temperature lapse rate provided\n');
  else
    fprintf(' - Calculate dew point temperature lapse rate by %s method\n',Parameter.LRTp_Vb);
  end
elseif strcmp(Forcing.IVTp_Vb,'Specific')
  fprintf('3. Dew point temperature calculated by Hum_Cal.m with specific humidity provided\n');
  fprintf(' - Calculate dew point temperature lapse rate by %s method\n',Parameter.LRTp_Vb);
elseif strcmp(Forcing.IVTp_Vb,'Relative')
  fprintf('3. Dew point temperature calculated by Hum_Cal.m with relative humidity provided\n');
  fprintf(' - Calculate dew point temperature lapse rate by %s method\n',Parameter.LRTp_Vb);
end

if strcmp(Parameter.Emis_DFn,'Off') % Emissivity provided or use built-in method
  fprintf('4. Calculate atmospheric emissivity by %s method\n',Parameter.EmTp_Vb);
else
  fprintf('4. Atmospheric emissivity provided\n')
end

if strcmp(Parameter.SMTp_Vb,'User') % Shadow mask provided or not
  fprintf('5. Shadow mask provided\n');
elseif strcmp(Parameter.SMTp_Vb,'Built-in')
  fprintf('5. Calculate shadow mask by CastShadow.m\n');
end

if strcmp(Parameter.SRTA_Fn,'Built-in') % Top-of-Atmoshpere shortwave provided or not
  fprintf('6. Calculate Top-of-Atmosphere shortwave by SW_DS.m line 34-38\n');
else
  fprintf('6. Top-of-Atmosphere shortwave provided\n');
end

if ~strcmp(Parameter.BSAl_DFn,'Off') && ~strcmp(Parameter.WSAl_DFn,'Off')
  fprintf('7. Black- and white-sky albedo provided\n'); % Black- and white-sky albedo provided or not
elseif strcmp(Parameter.BSAl_DFn,'Off') && strcmp(Parameter.WSAl_DFn,'Off')
  fprintf('7. Linearly interpolate input albedo at original resolution\n');
end

if strcmp(Forcing.Roug_Fn,'Off') % Surface roughness considered or not
  fprintf('8. Neglect surface roughness and displacement height\n');
else
  fid=fopen(Parameter.LCRg_Fn);
  LCRg1=fgets(fid);
  LCRg1=strsplit(LCRg1,',');
  LCRg=textscan(fid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %s','Delimiter',',');
  fclose(fid);
  LCRg=[LCRg1;LCRg];
  if strcmp(Forcing.DisH_Fn,'Off') % Displacement height considered or not
    fprintf('8. Neglect displacement height\n');
  else
    fprintf('8. Surface roughness and displacement height provided\n');
  end
  if strcmp(Parameter.VgIx_DFn,'Off') % Vegetation index provided or not
    fprintf(' - Neglect relative vegetation index\n');
    Vdn=[];
  else
    fprintf(' - Relative vegetation index provided\n');
    fid=fopen(Parameter.VITL_Fn);
    Vdn=textscan(fid,'%s','EndOfLine','\r\n');
    fclose(fid);
    Vdn=datenum(Vdn{1},Parameter.VIFm_Vb); % Datenum list of VI records
  end
end

if strcmp(Forcing.WsTp_Vb,'Component') % Component or total wind
  fprintf('9. U- and V-component of wind provided\n');
elseif strcmp(Forcing.WsTp_Vb,'Total')
  fprintf('9. Total wind provided\n');
end

%% Execute the calculation
fprintf('Execute the calculation\n');
stT=datenum(Domain.StTm_Vb,Domain.TmFm_Vb);
edT=datenum(Domain.EdTm_Vb,Domain.TmFm_Vb);
Trs=(edT-stT)/(Domain.NTmS_Vb-1);
T=stT:Trs:edT;
fprintf('Time resolution of forcing datasets is %s\n',duration(days(Trs),'Format','hh:mm:ss'));
clear LCRg1 fid stT edT Trs

if ~pflg
  for t=1:length(T)
    AtmDS_exe_sub(Domain,Forcing,Parameter,Output,T(t),LCRg,Vdn,flg);
  end

else
  parfor t=1:length(T)
    AtmDS_exe_sub(Domain,Forcing,Parameter,Output,T(t),LCRg,Vdn,flg);
  end
end
end

function AtmDS_exe_sub(Domain,Forcing,Parameter,Output,dn,LCRg,Vdn,flg)
Ds=datestr(dn,Domain.TmFm_Vb);
fprintf('Time step %s started',Ds);
[ys,~,~,~,~,~]=datevec(Ds,Domain.TmFm_Vb);
DoY=date2doy(dn);

if ischar(Forcing.MHti_Vb)
  MHti_b=retrFn(Forcing.MHti_Vb,Domain.TmFm_Vb,dn,ys);
  MHti=readCls(MHti_b);
  MHti=MHti(2:end-1,2:end-1);
else
  MHti_b=Forcing.MHti_Vb;
  MHti=MHti_b;
end

%% Downscaling of atmospheric variables
% Lapse rate for air temperature
Tair_bFn=retrFn(Forcing.Tair_bFn,Domain.TmFm_Vb,dn,ys);
if strcmp(Parameter.LRTp_Vb,'User')
  LRTa_Fn=retrFn(Parameter.LRTa_Fn,Parameter.LRFm_Vb,dn,ys);
  Tair=readCls(Tair_bFn);
  Tair=Tair(2:end-1,2:end-1);
  Tair=Temp_Adj(Tair,MHti,Forcing.MHto_Vb,LRTa_Fn);
else
  [LRTa,Tair,~]=LR_Temp(Parameter.LRTp_Vb,Tair_bFn,MHti_b,Domain.Elev_bFn,Forcing.MHto_Vb);
end
if ~strcmp(Output.Tair_OFn,'Off')
  fn=retrFn(Output.Tair_OFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Tair');
end

% Adjust height of pressure
Psfc_bFn=retrFn(Forcing.Psfc_bFn,Domain.TmFm_Vb,dn,ys);
if flg
  Psfc=Pair_Adj(Psfc_bFn,MHti_b,Tair_bFn,Forcing.MHto_Vb);
else
  Psfc=readCls(Psfc_bFn);
end
Psfc=Psfc(2:end-1,2:end-1);
if ~strcmp(Output.Psfc_OFn,'Off')
  fn=retrFn(Output.Psfc_OFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Psfc');
end

% Lapse rate for dew point temperature
if strcmp(Forcing.IVTp_Vb,'DewPoint')
  Tdew_bFn=retrFn(Forcing.IVar_bFn,Domain.TmFm_Vb,dn,ys);
  if strcmp(Parameter.LRTp_Vb,'User')
    LRTd_Fn=retrFn(Parameter.LRTd_Fn,Parameter.LRFm_Vb,dn,ys);
    Tdew=readCls(Tdew_bFn);
    Tdew=Tdew(2:end-1,2:end-1);
    Tdew=Temp_Adj(Tdew,MHti,Forcing.MHto_Vb,LRTd_Fn);
  else
    [LRTd,Tdew,~]=LR_Temp(Parameter.LRTp_Vb,Tdew_bFn,MHti_b,Domain.Elev_bFn,Forcing.MHto_Vb);
  end
else
  Hum_bFn=retrFn(Forcing.IVar_bFn,Domain.TmFm_Vb,dn,ys);
  Tdew_bFn=Hum_Cal(Tair_bFn,Psfc_bFn,Forcing.IVTp_Vb,Hum_bFn);
  [LRTd,Tdew,~]=LR_Temp(Parameter.LRTp_Vb,Tdew_bFn,MHti_b,Domain.Elev_bFn,Forcing.MHto_Vb);
end
if ~strcmp(Output.Tdew_OFn,'Off')
  fn=retrFn(Output.Tdew_OFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Tdew');
end

% Calculate humidity
[RHum,SHum]=Hum_Cal(Tair,Psfc,'DewPoint',Tdew);
if ~strcmp(Output.RHum_OFn,'Off')
  fn=retrFn(Output.RHum_OFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'RHum');
end
if ~strcmp(Output.SHum_OFn,'Off')
  fn=retrFn(Output.SHum_OFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'SHum');
end
fprintf('-');

% Atmospheric variables
Elev=readCls(Domain.Elev_bFn);
Elev=Elev(2:end-1,2:end-1);
[Tair_D,Tdew_D,Psfc_D,SHum_D,RHum_D]=AtmFrc_DS(Tair,LRTa,Tdew,LRTd,Psfc,Elev,Domain.Elev_DFn);
Tair_DFn=retrFn(Output.Tair_DFn,Domain.TmFm_Vb,dn,ys);
save(Tair_DFn,'Tair_D');
Tdew_DFn=retrFn(Output.Tdew_DFn,Domain.TmFm_Vb,dn,ys);
save(Tdew_DFn,'Tdew_D');
Psfc_DFn=retrFn(Output.Psfc_DFn,Domain.TmFm_Vb,dn,ys);
save(Psfc_DFn,'Psfc_D');
fn=retrFn(Output.SHum_DFn,Domain.TmFm_Vb,dn,ys);
save(fn,'SHum_D');
fn=retrFn(Output.RHum_DFn,Domain.TmFm_Vb,dn,ys);
save(fn,'RHum_D');
clear Tair_D Tdew_D Psfc_D SHum_D RHum_D LRTa LRTd RHum SHum Tair_bFn Psfc_bFn Hum_bFn Tdew_bFn MHti_b
clear Ds flg Elev
fprintf('-');

%% Downscaling of longwave radiation
LRad_Fn=retrFn(Forcing.LRad_Fn,Domain.TmFm_Vb,dn,ys);
if ~strcmp(Parameter.Emis_DFn,'Off')
  Emis_DFn=retrFn(Parameter.Emis_DFn,Domain.TmFm_Vb,dn,ys);
  LRad_D=LW_DS(LRad_Fn,Tair,Tdew,Tair_DFn,Tdew_DFn,Parameter.EmTp_Vb,Emis_DFn);
else
  LRad_D=LW_DS(LRad_Fn,Tair,Tdew,Tair_DFn,Tdew_DFn,Parameter.EmTp_Vb);
end
fn=retrFn(Output.LRad_DFn,Domain.TmFm_Vb,dn,ys);
save(fn,'LRad_D');
clear LRad_D LRad_Fn Emis_DFn Tdew Tair Tdew_DFn Tair_DFn
fprintf('-');

%% Downscaling of shortwave radiation
% Shadow
Azim_bFn=retrFn(Parameter.Azim_bFn,Domain.TmFm_Vb,dn,ys);
Azim_D=interp2d(Domain,Azim_bFn);
SoEl_bFn=retrFn(Parameter.SoEl_bFn,Domain.TmFm_Vb,dn,ys);
SoEl_D=interp2d(Domain,SoEl_bFn);
SdMk_DFn=retrFn(Parameter.SdMk_DFn,Domain.TmFm_Vb,dn,ys);
if strcmp(Parameter.SMTp_Vb,'User')
  if exist(SdMk_DFn,'file')==0
    [~,~,SdMk_D]=CastShadow(Domain.Xcor_DFn,Domain.Ycor_DFn,Domain.Lati_DFn,Domain.Elev_DFn,Azim_D,...
        SoEl_D);
  else
    SdMk_D=readCls(SdMk_DFn);
  end
elseif strcmp(Parameter.SMTp_Vb,'Built-in')
  if exist(SdMk_DFn,'file')==0
    [~,~,SdMk_D]=CastShadow(Domain.Xcor_DFn,Domain.Ycor_DFn,Domain.Lati_DFn,Domain.Elev_DFn,Azim_D,...
        SoEl_D);
    if ~isempty(find(SdMk_D==1, 1))
      save(SdMK_DFn,'SdMk');
    end
  else
    SdMk_D=readCls(SdMk_DFn);
  end
end
fprintf('-');

% Shortwave
SRad_Fn=retrFn(Forcing.SRad_Fn,Domain.TmFm_Vb,dn,ys);
if ~strcmp(Parameter.SRTA_Fn,'Built-in')
  Iflg='User';
  SRTA_Fn=retrFn(Parameter.SRTA_Fn,Domain.TmFm_Vb,dn,ys);
else
  Iflg='Built-in';
  SRTA_Fn=DoY;
end
Albd_Fn=retrFn(Forcing.Albd_Fn,Domain.TmFm_Vb,dn,ys);
if ~strcmp(Parameter.BSAl_DFn,'Off') && ~strcmp(Parameter.WSAl_DFn,'Off')
  BSAl_DFn=retrFn(Parameter.BSAl_DFn,Parameter.AlFm_Vb,dn,ys);
  WSAl_DFn=retrFn(Parameter.WSAl_DFn,Parameter.AlFm_Vb,dn,ys);
  SRad_D=SW_DS(SRad_Fn,Iflg,SRTA_Fn,Psfc,Psfc_DFn,Domain.Aspt_DFn,Domain.Slop_DFn,SdMk_D,...
      Domain.SVFa_DFn,Azim_D,SoEl_D,Albd_Fn,BSAl_DFn,WSAl_DFn);
elseif strcmp(Parameter.BSAl_DFn,'Off') && strcmp(Parameter.WSAl_DFn,'Off')
  SRad_D=SW_DS(SRad_Fn,Iflg,SRTA_Fn,Psfc,Psfc_DFn,Domain.Aspt_DFn,Domain.Slop_DFn,SdMk_D,...
      Domain.SVFa_DFn,Azim_D,SoEl_D,Albd_Fn);
end
fn=retrFn(Output.SRad_DFn,Domain.TmFm_Vb,dn,ys);
save(fn,'SRad_D');
clear SRad_D Azim_D SoEl_D SdMk_D Azim_bFn SoEl_bFn SdMk_DFn SRad_Fn SRTA_Fn Iflg Albd_Fn DoY BSAl_DFn
clear WSAl_DFn
fprintf('-');

%% Downscaling of wind speed
% Surface roughness and zero-plane displacement height
LCWt_DFn=retrFn(Parameter.LCWt_DFn,Parameter.LCFm_Vb,dn,ys);
if ~strcmp(Forcing.Roug_Fn,'Off')
  if ~strcmp(Parameter.VgIx_DFn,'Off')
    id=find(dn-Vdn>=0,1,'last');
    dnb=Vdn(id);
    fn=retrFn(Parameter.VgIx_DFn,Parameter.VIFm_Vb,dnb,ys); % VI before dn
    id=find(dn-Vdn<0,1,'first');
    dna=Vdn(id);
    fn1=retrFn(Parameter.VgIx_DFn,Parameter.VIFm_Vb,dna,ys); % VI after dn
    VgIx=Tinterp2D(fn,fn1,dnb,dna,dn);
  else
    VgIx=1;
  end
  Roug_Fn=retrFn(Forcing.Roug_Fn,Domain.TmFm_Vb,dn,ys);
  if ~strcmp(Forcing.DisH_Fn,'Off')
    DisH_Fn=retrFn(Forcing.DisH_Fn,Domain.TmFm_Vb,dn,ys);
    [Roug_D,DisH_D]=SfcRgh_DS(Roug_Fn,dn,Domain.Elev_DFn,LCWt_DFn,LCRg,VgIx,DisH_Fn);
  else
    DisH_Fn=0;
    Roug_D=SfcRgh_DS(Roug_Fn,dn,Domain.Elev_DFn,LCWt_DFn,LCRg,VgIx);
    DisH_D=0;
  end
else
  Roug_Fn=1;
  DisH_Fn=0;
  Roug_D=1;
  DisH_D=0;
end
fprintf('-');

% Adjust wind speed for friction velocity
if strcmp(Forcing.WsTp_Vb,'Component')
  Wspd_Fn=[retrFn(Forcing.Uspd_Fn,Domain.TmFm_Vb,dn,ys);retrFn(Forcing.Vspd_Fn,Domain.TmFm_Vb,dn,ys)];
  if ~strcmp(Output.Uspd_OFn,'Off') && ~strcmp(Output.Vspd_OFn,'Off')
    [~,~,Uspd,Vspd]=Wind_DS(Wspd_Fn,MHti,Forcing.MHto_Vb,Roug_Fn,Roug_Fn,DisH_Fn,DisH_Fn);
    fn=retrFn(Output.Uspd_OFn,Domain.TmFm_Vb,dn,ys);
    save(fn,'Uspd');
    fn=retrFn(Output.Vspd_OFn,Domain.TmFm_Vb,dn,ys);
    save(fn,'Vspd');
  end
  [~,~,Uspd_D,Vspd_D]=Wind_DS(Wspd_Fn,MHti,Forcing.MHto_Vb,Roug_Fn,Roug_D,DisH_Fn,DisH_D);
  fn=retrFn(Output.Uspd_DFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Uspd_D');
  fn=retrFn(Output.Vspd_DFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Vspd_D');

elseif strcmp(Forcing.WsTp_Vb,'Total')
  Wspd_Fn=retrFn(Forcing.Uspd_Fn,Domain.TmFm_Vb,dn,ys);
  if ~strcmp(Output.Uspd_OFn,'Off')
    Wspd=Wind_DS(Wspd_Fn,MHti,Forcing.MHto_Vb,Roug_Fn,Roug_Fn,DisH_Fn,DisH_Fn);
    fn=retrFn(Output.Uspd_OFn,Domain.TmFm_Vb,dn,ys);
    save(fn,'Wspd');
  end
  Wspd_D=Wind_DS(Wspd_Fn,MHti,Forcing.MHto_Vb,Roug_Fn,Roug_D,DisH_Fn,DisH_D);
  fn=retrFn(Output.Uspd_DFn,Domain.TmFm_Vb,dn,ys);
  save(fn,'Wspd_D');
end
fprintf('DONE\n');
end

function im_D=interp2d(Domain,fn)
Ycor=readCls(Domain.Ycor_bFn);
Xcor=readCls(Domain.Xcor_bFn);
im=readCls(fn);
Ycor_D=readCls(Domain.Ycor_DFn);
Xcor_D=readCls(Domain.Xcor_DFn);
im_D=interp2(Xcor,Ycor,im,Xcor_D,Ycor_D);
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

function fn=retrFn(fn,TsFm,dn,ys)
ds=datestr(dn,TsFm);
fn=strrep(fn,TsFm,ds);
if strcmp(fileparts(fn),'/yyyy') % Have year path or not
  fn=strrep(fn,'yyyy',ys);
end
end
