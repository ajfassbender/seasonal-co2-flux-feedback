%% Figure S7 Calculations
%
% DESCRIPTION:
%   Offline flux calculations using monthly model output.
%   Ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   flux_w92: CO2 flux calculation based on Wanninkhof (1992).
%
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;


%% Calculations -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

% Set file path
fpath = [data_path '/Ensemble_Members/m1/'];

co2     = reshape(ncread([fpath  'co2.nc'],'co2'),[],1);  
sst     = reshape(ncread([fpath  'sst.nc'],'sst'),[],1);  
sss     = reshape(ncread([fpath  'sss.nc'],'sss'),[],1);  
wind    = reshape(ncread([fpath  'wind.nc'],'wind'),[],1);  
co2_atm = reshape(ncread([fpath  'co2_atm.nc'],'co2_atm'),[],1);

g = find(wind >= 0 & sss > 0);

delta   = co2 - co2_atm;
[f,~,~] = flux_w92(wind(g),delta(g),sst(g),sss(g));
fx      = NaN(size(delta));
fx(g)   = f; clear f

flux_offline = reshape(fx,dim(1),dim(2),dim(3));% mol/m2/day

% Save variables to .nc files 
disp('Saving')  
nccreate([data_path '/Ensemble_Members/flux_offline_EM1.nc'],'flux_offline','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)})
ncwrite( [data_path '/Ensemble_Members/flux_offline_EM1.nc'],'flux_offline',flux_offline)
ncwriteatt([data_path '/Ensemble_Members/flux_offline_EM1.nc'],'flux_offline','units','mol/m2/d')
ncwriteatt([data_path '/Ensemble_Members/flux_offline_EM1.nc'],'flux_offline','long_name','Flux calculated offline from monthly co2, sst, sss, wind speed. onlu for Ensemble Member 1.')

%% Clean up

clear;close all;clc
