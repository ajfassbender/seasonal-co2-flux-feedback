%% Figure 8 Calculations 
%
% DESCRIPTION:
%   Calculate 2090s ensemble mean seasonal cycles.
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

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size;  

% Pre-Allocate
delta        = NaN(dim(1),dim(2),12,30);
delta_nl     = NaN(dim(1),dim(2),12,30);
winds        = NaN(dim(1),dim(2),12,30);
flux         = NaN(dim(1),dim(2),12,30);
flux_nl      = NaN(dim(1),dim(2),12,30);
delta_bp     = NaN(dim(1),dim(2),12,30);
delta_t      = NaN(dim(1),dim(2),12,30);
delta_bp_nl  = NaN(dim(1),dim(2),12,30);
delta_t_nl   = NaN(dim(1),dim(2),12,30);

st_0 = now;
for i = 1:30
    clearvars -except data_path i yri st_0 delta delta_nl winds flux flux_nl delta_bp delta_t delta_bp_nl delta_t_nl dim

    disp(['Start ensemble member ' num2str(i)])   
    st_1 = now; 
    
    fpath = [data_path '/Ensemble_Members/m' num2str(i)];
    
    co2_atm      = ncread([fpath  '/co2_atm.nc'],    'co2_atm',    [1 1 yri(1)], [Inf Inf length(yri)]);
    co2_corr     = ncread([fpath  '/co2_corr.nc'],   'co2_corr',   [1 1 yri(1)], [Inf Inf length(yri)]);
    co2_corr_nl  = ncread([fpath  '/co2_corr_nl.nc'],'co2_corr_nl',[1 1 yri(1)], [Inf Inf length(yri)]);
    wind         = ncread([fpath  '/wind.nc'],'wind', [1 1 yri(1)], [Inf Inf length(yri)]);
    flux_keff_nl = (365/12).*ncread([fpath  '/flux_keff_nl.nc'],'flux_keff_nl', [1 1 yri(1)], [Inf Inf length(yri)]);% mol/m2/mo
    flux_keff    = (365/12).*ncread([fpath  '/flux_keff.nc'],'flux_keff', [1 1 yri(1)], [Inf Inf length(yri)]);
 
    % Calculate bias-correced thermal and bp pCO2 terms
    co2_bias        = ncread([fpath  '/co2_bias.nc'], 'co2_bias', [1 1 yri(1)], [Inf Inf length(yri)]);
    co2_corr_bp     = ncread([fpath  '/co2_bp.nc'],   'co2_bp',   [1 1 yri(1)], [Inf Inf length(yri)]) + co2_bias;
    co2_corr_t      = ncread([fpath  '/co2_t.nc'],    'co2_t',    [1 1 yri(1)], [Inf Inf length(yri)]) + co2_bias;
    co2_corr_bp_nl  = ncread([fpath  '/co2_bp_nl.nc'],'co2_bp_nl',[1 1 yri(1)], [Inf Inf length(yri)]) + co2_bias;
    co2_corr_t_nl   = ncread([fpath  '/co2_t_nl.nc'], 'co2_t_nl', [1 1 yri(1)], [Inf Inf length(yri)]) + co2_bias;

    % 2090s Wind Speed, Flux, and Air-Sea pCO2 Disequilibrium Seasonal Cycles
    co2_atm_avg        = nanmean(reshape(co2_atm,         dim(1)*dim(2),12,10),3);
    delta(:,:,:,i)       = reshape(nanmean(reshape(co2_corr,dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    delta_nl(:,:,:,i)    = reshape(nanmean(reshape(co2_corr_nl,     dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    delta_bp(:,:,:,i)    = reshape(nanmean(reshape(co2_corr_bp,     dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    delta_t(:,:,:,i)     = reshape(nanmean(reshape(co2_corr_t,      dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    delta_bp_nl(:,:,:,i) = reshape(nanmean(reshape(co2_corr_bp_nl,  dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    delta_t_nl(:,:,:,i)  = reshape(nanmean(reshape(co2_corr_t_nl,   dim(1)*dim(2),12,10),3) - co2_atm_avg,dim(1),dim(2),12);
    winds(:,:,:,i)       = reshape(nanmean(reshape(wind,            dim(1)*dim(2),12,10),3),dim(1),dim(2),12); 
    flux(:,:,:,i)        = reshape(nanmean(reshape(flux_keff,       dim(1)*dim(2),12,10),3),dim(1),dim(2),12); 
    flux_nl(:,:,:,i)     = reshape(nanmean(reshape(flux_keff_nl,    dim(1)*dim(2),12,10),3),dim(1),dim(2),12); 
    
   disp(['Complete Ensemble Member ' num2str(i)])
   minutes((now - st_1)*24*60)
end
disp('Finish Calculations')
minutes((now - st_0)*24*60)
% ~12 minutes total


%% Compute Ensemble Mean Values -----------------------------------

delta_co2        = reshape(nanmean(delta,4),360,180,12);
delta_co2_nl     = reshape(nanmean(delta_nl,4),360,180,12);
delta_co2_bp     = reshape(nanmean(delta_bp,4),360,180,12);
delta_co2_t      = reshape(nanmean(delta_t,4),360,180,12);
delta_co2_bp_nl  = reshape(nanmean(delta_bp_nl,4),360,180,12);
delta_co2_t_nl   = reshape(nanmean(delta_t_nl,4),360,180,12);
wind_speed       = reshape(nanmean(winds,4),360,180,12);
co2_flux         = reshape(nanmean(flux,4),360,180,12);
co2_flux_nl      = reshape(nanmean(flux_nl,4),360,180,12);
     

%% Save variables to .nc files -----------------------------------
disp('Saving')

fpath = [data_path '/Ensemble_Members'];
[n,m,o] = size(co2_flux);

vars  = {'delta_co2','delta_co2_nl','wind_speed','co2_flux','co2_flux_nl','delta_co2_bp','delta_co2_t','delta_co2_bp_nl','delta_co2_t_nl'};
uni   = {'uatm','uatm','m/s','mol/m2/mo','mol/m2/mo','uatm','uatm','uatm','uatm'};
for i = 1:length(vars)
    uvar = char(vars(i));
    eval(['x = ' uvar ';'])
    nccreate([fpath   '/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'],uvar,'Dimensions',{'lon',n,'lat',m,'time',o})
    ncwrite([fpath    '/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'],uvar,x)
    ncwriteatt([fpath '/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'],uvar,'units',char(uni(i)))
    ncwriteatt([fpath '/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'],uvar,'long_name','2090s ensemble mean seasonal cycle')
end


%% -------- Clean Up

clear;close all;clc
