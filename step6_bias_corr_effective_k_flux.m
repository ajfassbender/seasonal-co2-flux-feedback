%% step6_bias_corr_effective_k_flux.m
%
% DESCRIPTION: 
%   Calculate effective gas transfer velocity (k_eff) from model output: (flux) / (delta pCO2).
%   Calculate smoothed decadal bias in reconstructed pCO2.
%   Offline calculation of sea-air CO2 flux using bias corrected pCO2 (co2_corr) and k_eff. 
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT: 
%    nonviable_keff (%) - ensemble mean percent of k_eff values globally that were not viable (~1.5%)
%    k_eff          (mol C/m2/d/uatm) - ensemble member effective gas transfer velocity
%    co2_bias       (microatmosphere) - ensemble member smoothed decadal bias in reconstructed pCO2
%    co2_corr       (microatmosphere) - ensemble member bias-corrected pCO2 values w/ chemical leverage
%    co2_corr_nl    (microatmosphere) - ensemble member bias-corrected pCO2 values w/o chemical leverage
%    co2_corr_tl    (microatmosphere) - ensemble member bias-corrected pCO2 values w/ thermal chemical leverage
%    co2_corr_bpl   (microatmosphere) - ensemble member bias-corrected pCO2 values w/ biophysical chemical leverage
%    flux_keff      (mol/m2/d) - ensemble member sea-air CO2 flux from co2_corr and k_eff
%    flux_keff_nl   (mol/m2/d) - ensemble member sea-air CO2 flux from co2_corr_nl and k_eff
%    flux_keff_tl   (mol/m2/d) - ensemble member sea-air CO2 flux from co2_corr_tl and k_eff
%    flux_keff_bpl  (mol/m2/d) - ensemble member sea-air CO2 flux from co2_corr_bpl and k_eff
%
% FUNCTIONS CALLED:
%   None. 
%
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your path')
data_path = cd;


%% Compute the pCO2 Bias, Effective Gas Transfer Velocity, and Flux -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

percent_nonviable_keff = NaN(30,1);

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim percent_nonviable_keff
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
        
    disp('Effective Gas Transfer Velocity Calculations')
    mflux = ncread([fpath  '/mflux.nc'],'mflux');   
    dco2  = ncread([fpath  '/dco2.nc'],'dco2');      
 	
    % Compute effective gas transfer velocity
	k_eff = -mflux./dco2; 
	clear mflux dco2 co2

	% Find outliers >10*median over all time
	med = median(k_eff,3);
	c   = find(abs(k_eff) > med*10);
	[c1 c2 c3] = ind2sub(size(k_eff),c);

	% Set outiers to the median value over all time
	for j = 1:length(c1)
		k_eff(c1(j),c2(j),c3(j)) = med(c1(j),c2(j));
	end	

	% Find remaining outliers >3*sigma from the mean over all time
	mn = abs(nanmean(k_eff,3)) + 3.*nanstd(k_eff,[],3);
	q  = find(abs(k_eff) > mn);
	[q1 q2 q3] = ind2sub(size(k_eff),q);

	% Set outiers to the mean value over all time
	for j = 1:length(q1)
		k_eff(q1(j),q2(j),q3(j)) = mn(q1(j),q2(j));
	end	

	% Determine the percentage of k_eff values that were not viable
	viable    = sum(sum(sum(isfinite(k_eff)))) - length(c1) - length(q1);
	nonviable = length(c1) + length(q1);
 	percent_nonviable_keff(ii) = nonviable./(viable+nonviable).*100; 
    clear viable nonviable q1 q2 q3 c1  c2 c3

    disp('CO2 Bias and Flux Calculations')
    co2      = ncread([fpath  '/co2.nc'],'co2');   
    co2_totl = ncread([fpath  '/co2_totl.nc'],'co2_totl'); 
    co2_nl   = ncread([fpath  '/co2_nl.nc'],'co2_nl'); 
    co2_atm  = ncread([fpath  '/co2_atm.nc'],'co2_atm'); 

    % Compute bias - smoothed over 10 years
    co2_bias = smoothdata(co2 - co2_totl,3,'lowess',120);  

	% Compute Fluxes
    co2_corr      = co2_totl + co2_bias;
    co2_corr_nl   = co2_nl   + co2_bias;
    flux_keff     = k_eff .* (co2_corr    - co2_atm);
    flux_keff_nl  = k_eff .* (co2_corr_nl - co2_atm);
    clear co2_atm co2_totl co2_nl

 	% Save variables to .nc files 
    disp('Saving')  

    nccreate([fpath   '/k_eff.nc'],'k_eff','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/k_eff.nc'],'k_eff',k_eff)
    ncwriteatt([fpath '/k_eff.nc'],'k_eff','units','molC/m2/d/uatm')
    ncwriteatt([fpath '/k_eff.nc'],'k_eff','long_name','effective gas transfer velocity')

    nccreate([fpath   '/co2_bias.nc'],'co2_bias','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_bias.nc'],'co2_bias',co2_bias)
    ncwriteatt([fpath '/co2_bias.nc'],'co2_bias','units','uatm')
    ncwriteatt([fpath '/co2_bias.nc'],'co2_bias','long_name','co2 bias (output - reconstruction)')

    nccreate([fpath   '/co2_corr.nc'],'co2_corr','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_corr.nc'],'co2_corr',co2_corr)
    ncwriteatt([fpath '/co2_corr.nc'],'co2_corr','units','uatm')
    ncwriteatt([fpath '/co2_corr.nc'],'co2_corr','long_name','reconstructed, bias corrected ocean co2 w/ chemical leverage')

    nccreate([fpath   '/co2_corr_nl.nc'],'co2_corr_nl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_corr_nl.nc'],'co2_corr_nl',co2_corr_nl)
    ncwriteatt([fpath '/co2_corr_nl.nc'],'co2_corr_nl','units','uatm')
    ncwriteatt([fpath '/co2_corr_nl.nc'],'co2_corr_nl','long_name','reconstructed, bias corrected ocean co2 w/o chemical leverage')

    nccreate([fpath   '/flux_keff.nc'],'flux_keff','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/flux_keff.nc'],'flux_keff',flux_keff)
    ncwriteatt([fpath '/flux_keff.nc'],'flux_keff','units','mol/m2/d')
    ncwriteatt([fpath '/flux_keff.nc'],'flux_keff','long_name','reconstructed, bias corrected co2 flux (w/ chemical leverage) using an effective gas transfer velocity')

    nccreate([fpath   '/flux_keff_nl.nc'],'flux_keff_nl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/flux_keff_nl.nc'],'flux_keff_nl',flux_keff_nl)
    ncwriteatt([fpath '/flux_keff_nl.nc'],'flux_keff_nl','units','mol/m2/d')
    ncwriteatt([fpath '/flux_keff_nl.nc'],'flux_keff_nl','long_name','reconstructed, bias corrected co2 flux (w/o chemical leverage) using an effective gas transfer velocity')
    
    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~6 min per ensemble member

nonviable_keff = nanmean(percent_nonviable_keff);   % 1.5%
nonviable_keff_stdev = std(percent_nonviable_keff); % 0.01%
nccreate([data_path   '/Ensemble_Members/nonviable_keff.nc'],'nonviable_keff')
ncwrite( [data_path   '/Ensemble_Members/nonviable_keff.nc'],'nonviable_keff',nonviable_keff)
ncwriteatt([data_path '/Ensemble_Members/nonviable_keff.nc'],'nonviable_keff','units','%')
ncwriteatt([data_path '/Ensemble_Members/nonviable_keff.nc'],'nonviable_keff','long_name','ensemble mean percent of effective gas transfer velocities that were not viable')


%% Compute fluxes for t and bp leverage -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
        
    k_eff    = ncread([fpath  '/k_eff.nc'],'k_eff');   
    co2_bias = ncread([fpath  '/co2_bias.nc'],'co2_bias'); 
    co2_atm  = ncread([fpath  '/co2_atm.nc'],'co2_atm');  
    co2_tl   = ncread([fpath  '/co2_tl.nc'],'co2_tl');   
    co2_bpl  = ncread([fpath  '/co2_bpl.nc'],'co2_bpl');     

	% pCO2 Bias Correction
    co2_corr_tl    = co2_tl  + co2_bias;	
    co2_corr_bpl   = co2_bpl + co2_bias;	
    clear co2_bias
    
    % Compute Fluxes % molC/m2/d
    flux_keff_tl    = k_eff .* (co2_corr_tl  - co2_atm);
    flux_keff_bpl   = k_eff .* (co2_corr_bpl - co2_atm);
    
    clear co2_atm 

 	% Save variables to .nc files 
    disp('Saving')  

    nccreate([fpath   '/co2_corr_tl.nc'],'co2_corr_tl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_corr_tl.nc'],'co2_corr_tl',co2_corr_tl)
    ncwriteatt([fpath '/co2_corr_tl.nc'],'co2_corr_tl','units','uatm')
    ncwriteatt([fpath '/co2_corr_tl.nc'],'co2_corr_tl','long_name','reconstructed, bias corrected ocean co2 w/ thermal leverage only')
   
    nccreate([fpath   '/co2_corr_bpl.nc'],'co2_corr_bpl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_corr_bpl.nc'],'co2_corr_bpl',co2_corr_bpl)
    ncwriteatt([fpath '/co2_corr_bpl.nc'],'co2_corr_bpl','units','uatm')
    ncwriteatt([fpath '/co2_corr_bpl.nc'],'co2_corr_bpl','long_name','reconstructed, bias corrected ocean co2 w/ biophysical leverage only')

    nccreate([fpath   '/flux_keff_tl.nc'],'flux_keff_tl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/flux_keff_tl.nc'],'flux_keff_tl',flux_keff_tl)
    ncwriteatt([fpath '/flux_keff_tl.nc'],'flux_keff_tl','units','mol/m2/d')
    ncwriteatt([fpath '/flux_keff_tl.nc'],'flux_keff_tl','long_name','reconstructed, bias corrected co2 flux (w/ thermal leverage only) using an effective gas transfer velocity')

    nccreate([fpath   '/flux_keff_bpl.nc'],'flux_keff_bpl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/flux_keff_bpl.nc'],'flux_keff_bpl',flux_keff_bpl)
    ncwriteatt([fpath '/flux_keff_bpl.nc'],'flux_keff_bpl','units','mol/m2/d')
    ncwriteatt([fpath '/flux_keff_bpl.nc'],'flux_keff_bpl','long_name','reconstructed, bias corrected co2 flux (w/ biophysical leverage only) using an effective gas transfer velocity')
    
    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~4 min per ensemble member


%% Clean up -----------------------------------

clear;close all;clc

