%% step9_delta_c_update_mean_pco2.m
%
% DESCRIPTION:
%   Calculate cumulative ocean carbon uptake using reconstructed and bias corrected annual mean pCO2 values.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT:    
%   EnsembleMean_C_Uptake (PgC) - file containing ensemble mean global, area-weighted and summed cumualtive ocean 
%                                 carbon uptake for fluxes based on annual mean pCO2 values:
%                                 'tot_uptake_annual_mean_co2', 'tot_uptake_annual_mean_co2_stdev'
%   EnsembleMean_C_Uptake (mol C/mo) - file containing ensemble mean global cumualtive ocean carbon uptake maps and standard
%                                      deviation maps for fluxes based on annual mean pCO2 values:
%                                      'uptake_annual_mean_co2', 'uptake_annual_mean_co2_stdev'
%
% FUNCTIONS CALLED:
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your path')
data_path = cd;


%% Calculations -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

area  = transpose(getArea(-89.5:1:89.5,.5:359.5));

% Pre-allocate
F_avg = NaN(dim(1),dim(2),30);
fluf_keff_avg_tot = NaN(30,1);

st_0 = now;
for ii = 1:30    
    clear x i_f f_keff_avg fluf_keff_avg_mo

    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];

    % Load Data (fluxes: mol/m2/d)
    i_f        = 1 - reshape(ncread( [fpath  '/ice_frac.nc'],'ice_frac'),[],1);  
    f_keff_avg = reshape(ncread( [fpath  '/flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2'),[],1);  
    
    % Apply ice mask
    x  = find(i_f > 0);
    f_keff_avg(x)  = i_f(x) .* f_keff_avg(x);      
    
    % Monthly flux (mol/m2/d --> mol/m2/mo)
    fluf_keff_avg_mo = reshape(f_keff_avg,360,180,1812) .* (365/12);

    % Cumulative flux over all time (mol/m2)
    F_avg(:,:,ii) = nansum(fluf_keff_avg_mo,3);

    % Area Weighted Global Cumulative Uptake (mol/m2 --> Pg C)
    fluf_keff_avg_tot(ii,1) = nansum(nansum(F_avg(:,:,ii) .* area)) .* 12 ./ (10^15);

    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)

%% Save Ensemble Mean -----------------------------------

fpath = [data_path '/Ensemble_Members'];

% Ensemble mean cumualtive carbon uptake % mol/m2 
uptake_annual_mean_co2         = nanmean(F_avg,3);
uptake_annual_mean_co2_stdev   = nanstd(F_avg,[],3);

% Ensemble mean area-weighted & summed cumulative carbon uptake (Pg C)
tot_uptake_annual_mean_co2       = nanmean(fluf_keff_avg_tot); 
tot_uptake_annual_mean_co2_stdev = nanstd(fluf_keff_avg_tot);

nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2','Dimensions',{'lon',dim(1),'lat',dim(2)})
ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2',uptake_annual_mean_co2)
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2','units','mol/m2')
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2','long_name','Global 1950-2100 cumulative ocean carbon uptake maps: using annual mean reconstructed & bias corrected co2')

nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2_stdev','Dimensions',{'lon',dim(1),'lat',dim(2)})
ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2_stdev',uptake_annual_mean_co2_stdev)
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2_stdev','units','mol/m2')
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2_stdev','long_name','Global 1950-2100 cumulative ocean carbon uptake stdev maps: using annual mean reconstructed & bias corrected co2')

nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2')
ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2',tot_uptake_annual_mean_co2)
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2','units','Pg C')
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2','long_name','Global, area-weighted and summed 1950-2100 cumulative ocean carbon uptake: using annual mean reconstructed & bias corrected co2')

nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2_stdev')
ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2_stdev',tot_uptake_annual_mean_co2_stdev)
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2_stdev','units','Pg C')
ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2_stdev','long_name','Global, area-weighted and summed 1950-2100 cumulative ocean carbon uptake stdev: using annual mean reconstructed & bias corrected co2')


%% Clean up -----------------------------------

clear;close all;clc
