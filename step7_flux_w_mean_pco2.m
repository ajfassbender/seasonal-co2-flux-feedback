%% step7_flux_w_mean_pco2.m
%
% DESCRIPTION:
%   Calculate sea-air CO2 fluxes using annual mean pCO2 values and monthly effective gas transfer velocities. 
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT: 
%    flux_keff_annual_mean_co2 (mol/m2/d) - ensemble member sea-air CO2 flux from annual mean pCO2 values and k_eff
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


%% Compute Flux using Annual Mean pCO2 Values that Have Been Reconstructed and Bias Corrected -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
        
    co2_atm  = ncread([fpath  '/co2_atm.nc'],'co2_atm');  
    co2      = ncread([fpath  '/co2.nc'],'co2');  
    k_eff    = ncread([fpath  '/k_eff.nc'],'k_eff');  
    
    % Calculate annual mean values, repeat annual value over 12 months of each year 
    co2_av = reshape(repmat(nanmean(reshape(co2,dim(1),dim(2),12,dim(3)/12),3),1,1,12,1),dim(1),dim(2),dim(3));
    flux_keff_annual_mean_co2 = k_eff .* (co2_av - co2_atm);
    
    % Save variables to .nc files 
    disp('Saving')  

    nccreate([fpath   '/flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2',flux_keff_annual_mean_co2)
    ncwriteatt([fpath '/flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2','units','mol/m2/d')
    ncwriteatt([fpath '/flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2','long_name','co2 flux (using annual mean co2 vlues w/ bias correction) recomputed using an effective gas transfer velocity')

    disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Clean up -----------------------------------

%clear;close all;clc
