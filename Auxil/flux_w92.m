function [Flux_estimate,w_k_pv,w_Co] = flux_w92(wind_speed,delta_pco2,SST,SSS)
%%
% [Flux_estimate] = flux_w92(wind_speed,delta_pco2,SST,SSS)
%
% units = mol/m2 day
%
% Wanninkhof, Rik (1992), Relationship between wind speed and gas exchange 
% over the ocean, JGR, 97, doi:10.1029/92JC00188.
%
% Input:
%   wind speed (m/s)
%   dpCO2 (uatm) - sea-air
%   sst (°C)	
%   sss (none)
  
% Schmidt # (Table A1 Wanninkhof, 1992)
aa = 2073.1;
bb = -125.62;
cc = 3.6276;
dd = -0.043219;
Sch = aa + bb.*SST + cc.*(SST.^2) + dd.*(SST.^3); % unitless

% Wind speed squared ((m/s)^2)
u2 = wind_speed.^2; 

% Gas transfer coefficient (cm/hr)
k = 0.31.*u2.*(Sch./660).^-0.5; 

% CO2 solubility Wanninkhof, 2014 --> from Weis 1974 (T in Kelvin)
% NOTE: Wanninkhof 1993 does provide a different version of this where the Ko value is in mol/kg/atm
T  = SST + 273.15; % C to kelvin
A1 = -58.0931; 
A2 =  90.5069;
A3 =  22.2940;
B1 =   0.027766;
B2 =  -0.025888;
B3 =   0.0050578;

lnKo = A1 + A2.*(100./T) + A3.*log(T./100)+ SSS.*(B1 + B2.*(T./100) + B3.*(T./100).^2); % mol/L/atm
Ko   = exp(lnKo);
Co   = Ko .* 1000; % mol/m3 atm 

Flux_estimate = Co.*k.*(10^-6.*delta_pco2).*(24/100);%mol/m3 atm * cm/hr * atm = mol/m3 *cm/hr * 1m/100cm * 24hr/d = mol/m2 day
w_k_pv = k;
w_Co   = Co;

end
