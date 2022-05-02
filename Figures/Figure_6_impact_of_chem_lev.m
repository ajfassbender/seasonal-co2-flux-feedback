%% Figure 6 Calculations
%
% DESCRIPTION:
%   Impact of chemical leverage on seasonal pCO2 values and fluxes.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   None.
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;


%% Make caluclations -----------------------------------

fpath = [data_path '/Ensemble_Members'];
lat   = ncread([fpath '/lat.nc'],  'lat');   
lon   = ncread([fpath '/lon.nc'],  'lon'); 
time  = ncread([fpath '/time.nc'], 'time'); 

yr   = str2num(datestr(time,'yyyy'));
yri  = find(2090 <= yr & yr < 2100);
uyr  = unique(yr);
uyri = find(2090 <= uyr & uyr <2100);

% Sesaonal indices for hemispheres
nh_s = 5:10;
nh_w = [1:4 11:12];

flux_s  = NaN(360,180,30);
flux_w  = NaN(360,180,30);
dco2_s  = NaN(360,180,30);
dco2_w  = NaN(360,180,30);
daflux  = NaN(360,180,30);

st_0 = now;
for i = 1:30
   clearvars -except data_path yr yri uyr uyri mo i st_0 lon lat flux_s flux_w dco2_s dco2_w daflux nh_s nh_w

   disp(['Start ensemble member ' num2str(i)])   
   st_1 = now; 
    
   fpath = [data_path '/Ensemble_Members/m' num2str(i)];
      
   % Load variables (and convert daily fluxes to monthly)
   flux_nl     = (365/12).*ncread([fpath  '/flux_keff_nl.nc'],'flux_keff_nl', [1 1 yri(1)], [Inf Inf length(yri)]); % mol/m2/mo
   flux        = (365/12).*ncread([fpath  '/flux_keff.nc'],'flux_keff', [1 1 yri(1)], [Inf Inf length(yri)]); % mol/m2/mo
   co2_corr_nl = ncread([fpath  '/co2_corr_nl.nc'],'co2_corr_nl', [1 1 yri(1)], [Inf Inf length(yri)]); 
   co2_corr    = ncread([fpath  '/co2_corr.nc'],'co2_corr', [1 1 yri(1)], [Inf Inf length(yri)]);

   % 2090s averages
   flux_nhs = NaN(360,180,length(uyri));
   flux_nhw = NaN(360,180,length(uyri));
   dco2_nhs = NaN(360,180,length(uyri));
   dco2_nhw = NaN(360,180,length(uyri));
   aflux    = NaN(360,180,length(uyri));

   % Cycle through each year
   for j = 1:length(uyri)
      clear q;q  = find(yr(yri) == uyr(uyri(j)));

      % Mean seasonal difference betwteen chem-leverage and no-leverage pCO2
      clear X; X       = co2_corr(:,:,q) - co2_corr_nl(:,:,q);
      dco2_nhs(:,:,j)  = nanmean(X(:,:,nh_s),3);
      dco2_nhw(:,:,j)  = nanmean(X(:,:,nh_w),3);

      % Annual uptake difference betwteen chem-leverage and no-leverage pCO2
      aflux(:,:,j)     = nansum(flux(:,:,q) - flux_nl(:,:,q),3);% mol/m2/mo --> mol/m2/yr
      
      % Seasonal cumulative flux difference betwteen chem-leverage and no-leverage 
      clear X; X       = flux(:,:,q) - flux_nl(:,:,q); % mol/m2/mo 
      flux_nhs(:,:,j)  = nansum(X(:,:,nh_s),3); % mol/m2/season 
      flux_nhw(:,:,j)  = nansum(X(:,:,nh_w),3); % mol/m2/season 
   end
   
   % Combine hemispheres and average over 2090s
   qn = find(lat>0);  
   qs = find(lat<0);
   flux_w(:,:,i) = nanmean([flux_nhs(:,qs,:) flux_nhw(:,qn,:)],3); % mol/m2/season  
   flux_s(:,:,i) = nanmean([flux_nhw(:,qs,:) flux_nhs(:,qn,:)],3); % mol/m2/season  
   dco2_w(:,:,i) = nanmean([dco2_nhs(:,qs,:) dco2_nhw(:,qn,:)],3); % uatm
   dco2_s(:,:,i) = nanmean([dco2_nhw(:,qs,:) dco2_nhs(:,qn,:)],3); % uatm
   daflux(:,:,i) = nanmean(aflux,3); % mol/m2/yr
   
   disp(['Complete Ensemble Member ' num2str(i)])
   minutes((now - st_1)*24*60)
end
disp('Finish Calculations')
minutes((now - st_0)*24*60)


%% Compute Ensemble Mean Values -----------------------------------

diff_co2_summer   = nanmean(dco2_s,3); % sesaonal mean pCO2 difference
diff_co2_winter   = nanmean(dco2_w,3);
diff_flux_summer  = nanmean(flux_s,3); % seasonal cumulative flux difference
diff_flux_winter  = nanmean(flux_w,3);
diff_aflux        = nanmean(daflux,3); % annual cumulative flux difference


%% Save variables to .nc files -----------------------------------
disp('Saving') 

fpath = [data_path '/Ensemble_Members'];
[n,m] = size(diff_aflux);

vars  = {'diff_co2_summer','diff_co2_winter','diff_flux_summer','diff_flux_winter','diff_aflux'};
uni   = {'uatm','uatm','mol/m2/season','mol/m2/season','mol/m2/year'};
for i = 1:length(vars)
    uvar = char(vars(i));
    eval(['x = ' uvar ';'])
    nccreate([fpath   '/Figure6_impact_of_chem_lev.nc'],uvar,'Dimensions',{'lon',n,'lat',m})
    ncwrite([fpath    '/Figure6_impact_of_chem_lev.nc'],uvar,x)
    ncwriteatt([fpath '/Figure6_impact_of_chem_lev.nc'],uvar,'units',char(uni(i)))
    ncwriteatt([fpath '/Figure6_impact_of_chem_lev.nc'],uvar,'long_name','ensemble mean difference between chem-leverage and no-leverage variable')
end


%% Clean Up -----------------------------------

clear;close all;clc
