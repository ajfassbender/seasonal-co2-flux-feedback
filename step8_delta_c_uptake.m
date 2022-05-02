%% step8_delta_c_uptake.
%
% DESCRIPTION:
%   Calculate cumulative ocean carbon uptake using reconstructed and bias corrected pCO2 values for
%   each chemical leverage configuration.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT:    
%   Cumulative_C_Uptake (PgC) - file containing global, area-weighted and summed cumualtive ocean carbon uptake 
%                               for model output (mflux) and each chemical leverage configuration: 'mflux_tot', 
%                               'flux_keff_tot', 'flux_keff_nl_tot', 'flux_keff_tl_tot', 'flux_keff_bpl_tot'
%   EnsembleMean_C_Uptake (PgC) - file containing ensemble mean global, area-weighted and summed cumualtive ocean 
%                                 carbon uptake for model output (mflux) and each chemical leverage configuration:
%                                 'tot_uptake_nl', 'tot_uptake_l', 'tot_uptake_tl', 'tot_uptake_bpl', 'tot_uptake_model'
%   EnsembleMean_C_Uptake (mol C/mo) - file containing ensemble mean global cumualtive ocean carbon uptake maps and standard
%                                      deviation maps for model output (mflux) and each chemical leverage configuration:
%                                      'uptake_nl', 'uptake_l', 'uptake_tl','uptake_bpl', 'uptake_model'
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


%% Calculate Carbon Uptake  -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 
area  = transpose(getArea(-89.5:1:89.5,.5:359.5));

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim  area
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    
    % Load Data (fluxes: mol/m2/d)
    i_f       = 1 - reshape(ncread( [fpath  '/ice_frac.nc'],'ice_frac'),[],1);   
    f_mflux   = reshape(-ncread([fpath  '/mflux.nc'],'mflux'),[],1);
    f_keff    = reshape(ncread( [fpath  '/flux_keff.nc'],'flux_keff'),[],1); 
    f_keff_nl = reshape(ncread( [fpath  '/flux_keff_nl.nc'],'flux_keff_nl'),[],1);

    % Apply ice mask
    x = find(i_f > 0);
    f_mflux(x)   = i_f(x) .* f_mflux(x); 
    f_keff(x)    = i_f(x) .* f_keff(x);
    f_keff_nl(x) = i_f(x) .* f_keff_nl(x);     
    
    % Monthly flux (mol/m2/d --> mol/m2/mo --> mol/m2)
    mflux_mo         = nansum(reshape(f_mflux,  dim(1),dim(2),dim(3)) .* (365/12),3);
    flux_keff_mo     = nansum(reshape(f_keff,   dim(1),dim(2),dim(3)) .* (365/12),3);
    flux_keff_nl_mo  = nansum(reshape(f_keff_nl,dim(1),dim(2),dim(3)) .* (365/12),3);

    % Area Weighted Global Cumulative Uptake (mol/m2 --> Pg C)
    mflux_tot        = nansum(nansum(mflux_mo        .* area)) .* 12 ./ (10^15);
    flux_keff_tot    = nansum(nansum(flux_keff_mo    .* area)) .* 12 ./ (10^15);
    flux_keff_nl_tot = nansum(nansum(flux_keff_nl_mo .* area)) .* 12 ./ (10^15);
    
    % Save variables to .nc files 
    disp('Saving')  

    nccreate([fpath   '/Cumulative_C_Uptake.nc'],'mflux_tot')
    ncwrite( [fpath   '/Cumulative_C_Uptake.nc'],'mflux_tot',mflux_tot)
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'mflux_tot','units','Pg C')
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'mflux_tot','long_name','Area-weighted 1950-2100 cumulative ocean carbon uptake from model output flux')

    nccreate([fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_tot')
    ncwrite( [fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_tot',flux_keff_tot)
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_tot','units','Pg C')
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_tot','long_name','Area-weighted 1950-2100 cumulative ocean carbon uptake from recomputed fluxes using monthly effective gas transfer velocities')

    nccreate([fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_nl_tot')
    ncwrite( [fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_nl_tot',flux_keff_nl_tot)
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_nl_tot','units','Pg C')
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_nl_tot','long_name','Area-weighted 1950-2100 cumulative ocean carbon uptake from recomputed fluxes (no chemical leverage) using monthly effective gas transfer velocities')

    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Add in the t and bp leverage values -----------------------------------

clearvars -except data_path 
clc

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 
area  = transpose(getArea(-89.5:1:89.5,.5:359.5));

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim area
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    
    % Load Data (mol/m2/d)
    i_f        = 1 - reshape(ncread( [fpath  '/ice_frac.nc'],'ice_frac'),[],1);   
    f_keff_bpl = reshape(ncread( [fpath  '/flux_keff_bpl.nc'],'flux_keff_bpl'),[],1); 
    f_keff_tl  = reshape(ncread( [fpath  '/flux_keff_tl.nc'],'flux_keff_tl'),[],1);

    % Apply ice mask
    x = find(i_f > 0);
    f_keff_bpl(x) = i_f(x) .* f_keff_bpl(x);  
    f_keff_tl(x)  = i_f(x) .* f_keff_tl(x);     
    
    % Monthly flux (mol/m2/d --> mol/m2/mo --> mol/m2)
    flux_keff_tl_mo  = nansum(reshape(f_keff_tl, dim(1),dim(2),dim(3)) .* (365/12),3);
    flux_keff_bpl_mo = nansum(reshape(f_keff_bpl,dim(1),dim(2),dim(3)) .* (365/12),3);

    % Area Weighted Global Cumulative Uptake (mol/m2 --> Pg C)
    flux_keff_tl_tot   = nansum(nansum(flux_keff_tl_mo  .* area)) .* 12 ./ (10^15);
    flux_keff_bpl_tot  = nansum(nansum(flux_keff_bpl_mo .* area)) .* 12 ./ (10^15);

    % Save variables to .nc files 
    disp('Saving')  

    nccreate([fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_tl_tot')
    ncwrite( [fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_tl_tot',flux_keff_tl_tot)
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_tl_tot','units','Pg C')
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_tl_tot','long_name','Area-weighted 1950-2100 cumulative ocean carbon uptake from recomputed fluxes (thermal leverage only) using monthly effective gas transfer velocities')

    nccreate([fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_bpl_tot')
    ncwrite( [fpath   '/Cumulative_C_Uptake.nc'],'flux_keff_bpl_tot',flux_keff_bpl_tot)
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_bpl_tot','units','Pg C')
    ncwriteatt([fpath '/Cumulative_C_Uptake.nc'],'flux_keff_bpl_tot','long_name','Area-weighted 1950-2100 cumulative ocean carbon uptake from recomputed fluxes (biophysical leverage only) using monthly effective gas transfer velocities')

    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Save Ensemble Mean Result -----------------------------------

clearvars -except data_path 
clc

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

mflux_map    = NaN(dim(1),dim(2),30);
flux_tl_map  = NaN(dim(1),dim(2),30);
flux_bpl_map = NaN(dim(1),dim(2),30);
flux_map     = NaN(dim(1),dim(2),30);
flux_nl_map  = NaN(dim(1),dim(2),30);    

st_0 = now;
for ii = 1:30    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];

    % Cumulative carbon uptake (mol/m2/d --> mol/m2/mo --> mol/m2)
    mflux_map(:,:,ii)     = sum((365/12) .* ncread([fpath '/mflux.nc'],'mflux'),3); 
    flux_tl_map(:,:,ii)   = sum((365/12) .* ncread([fpath '/flux_keff_tl.nc'], 'flux_keff_tl'),3); 
    flux_bpl_map(:,:,ii)  = sum((365/12) .* ncread([fpath '/flux_keff_bpl.nc'],'flux_keff_bpl'),3);
    flux_map(:,:,ii)      = sum((365/12) .* ncread([fpath '/flux_keff.nc'],    'flux_keff'),3);    
    flux_nl_map(:,:,ii)   = sum((365/12) .* ncread([fpath '/flux_keff_nl.nc'], 'flux_keff_nl'),3);

    % Area-weighted and summed total cumulative carbon uptake (Pg C)
    mflux_tot(ii,1)         = ncread([fpath  '/Cumulative_C_Uptake.nc'],'mflux_tot');         % model output flux
    flux_keff_tl_tot(ii,1)  = ncread([fpath  '/Cumulative_C_Uptake.nc'],'flux_keff_tl_tot');  % thermal leverage only reconstructed co2
    flux_keff_bpl_tot(ii,1) = ncread([fpath  '/Cumulative_C_Uptake.nc'],'flux_keff_bpl_tot'); % bp leverage only reconstructed co2
    flux_keff_nl_tot(ii,1)  = ncread([fpath  '/Cumulative_C_Uptake.nc'],'flux_keff_nl_tot');  % no leverage reconstructed co2
    flux_keff_tot(ii,1)     = ncread([fpath  '/Cumulative_C_Uptake.nc'],'flux_keff_tot');     % total chem leverage reconstructed co2

    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)

% Ensemble mean cumualtive carbon uptake % mol/m2 
uptake_nl        = nanmean(flux_nl_map,3); 
uptake_nl_sd     = nanstd(flux_nl_map,[],3);
uptake_l         = nanmean(flux_map,3);
uptake_l_sd      = nanstd(flux_map,[],3);
uptake_tl        = nanmean(flux_tl_map,3);
uptake_tl_sd     = nanstd(flux_tl_map,[],3);
uptake_bpl       = nanmean(flux_bpl_map,3);
uptake_bpl_sd    = nanstd(flux_bpl_map,[],3);
uptake_model     = nanmean(mflux_map,3);
uptake_model_sd  = nanstd(mflux_map,[],3);

% Ensemble mean area-weighted & summed cumulative carbon uptake (Pg C)
tot_uptake_nl        = nanmean(flux_keff_nl_tot);
tot_uptake_nl_sd     = nanstd(flux_keff_nl_tot);
tot_uptake_l         = nanmean(flux_keff_tot);
tot_uptake_l_sd      = nanstd(flux_keff_tot);
tot_uptake_tl        = nanmean(flux_keff_tl_tot);
tot_uptake_tl_sd     = nanstd(flux_keff_tl_tot);
tot_uptake_bpl       = nanmean(flux_keff_bpl_tot);
tot_uptake_bpl_sd    = nanstd(flux_keff_bpl_tot);
tot_uptake_model     = nanmean(mflux_tot);
tot_uptake_model_sd  = nanstd(mflux_tot);

% Save Ensemble Mean 
fpath = [data_path '/Ensemble_Members'];

vname = {'uptake_nl','uptake_l','uptake_tl','uptake_bpl','uptake_model'};
lname = {'no leverage','leverage','thermal leverage','biophysical leverage','model output'};

for i = 1:5
    var = char(vname(i));  
    var_stdev = [var '_sd'];

    nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],var,'Dimensions',{'lon',dim(1),'lat',dim(2)})
    ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],var, eval(var))
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],var,'units','mol/m2')
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],var,'long_name',['Global 1950-2100 cumulative ocean carbon uptake maps: ' char(lname(i))])
    
    nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],var_stdev,'Dimensions',{'lon',dim(1),'lat',dim(2)})
    ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],var_stdev,eval(var_stdev))
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],var_stdev,'units','mol/m2')
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],var_stdev,'long_name',['Global 1950-2100 cumulative ocean carbon uptake stdev maps: ' char(lname(i))])
    
    nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],['tot_' var])
    ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],['tot_' var],eval(['tot_' var]))
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],['tot_' var],'units','Pg C')
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],['tot_' var],'long_name',['Global, area-weighted and summed 1950-2100 cumulative ocean carbon uptake: ' char(lname(i))])

    nccreate([fpath   '/EnsembleMean_C_Uptake.nc'],['tot_' var_stdev])
    ncwrite( [fpath   '/EnsembleMean_C_Uptake.nc'],['tot_' var_stdev],eval(['tot_' var_stdev]))
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],['tot_' var_stdev],'units','Pg C')
    ncwriteatt([fpath '/EnsembleMean_C_Uptake.nc'],['tot_' var_stdev],'long_name',['Global, area-weighted and summed 1950-2100 cumulative ocean carbon uptake stdev: ' char(lname(i))])
end

% Ensemble Mean Uptake in Pg C 

100 .* (tot_uptake_l - tot_uptake_nl)  ./ tot_uptake_l % globally  8% greater transient uptake
100 .* (tot_uptake_l - tot_uptake_tl)  ./ tot_uptake_l % globally  9% greater tlev uptake
100 .* (tot_uptake_l - tot_uptake_bpl) ./ tot_uptake_l % globally 17% greater transient uptake

%% Clean up -----------------------------------

clear;close all;clc
