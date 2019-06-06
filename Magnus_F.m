function E=Magnus_F(T)
abs0=-273.15;

% Coefficient for Magnus formula adpoted from Buck (1981)
Aw=611.21; % If Ta>5, use coeffcient of the ew2 curve in Buck (1981)
Bw=17.368;
Cw=238.88;
Am=611.21; % If -5<=Ta<=5, use coeffcient of the ew1 curve in Buck (1981)
Bm=17.502;
Cm=240.97;
Ai=611.15; % If Ta<-5, use coeffcient of the ei2 curve in Buck (1981)
Bi=22.452;
Ci=272.55;

% Saturated vapor pressure (Pa)
E=Aw*exp(Bw*(T+abs0)./(T+abs0+Cw)); % Magnus formula
E(T<=5-abs0)=Am*exp(Bm*(T(T<=5-abs0)+abs0)./(T(T<=5-abs0)+abs0+Cm));
E(T<-5-abs0)=Ai*exp(Bi*(T(T<-5-abs0)+abs0)./(T(T<-5-abs0)+abs0+Ci));
end
