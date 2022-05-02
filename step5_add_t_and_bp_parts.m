%% step5_add_t_and_bp_parts.m
% 
% DESCRIPTION:
%   Sum thermal and biophysical pCO2 components.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT:    
%   co2_nl   (microatmosphere) - ensemble member reconstructed pCO2 w/o chemical leverage (nl)
%   co2_bpl  (microatmosphere) - ensemble member reconstructed pCO2 w/ biophysical chemical leverage (bpl)
%   co2_tl   (microatmosphere) - ensemble member reconstructed pCO2 w/ thermal chemical leverage (tl)
%   co2_totl (microatmosphere) - ensemble member reconstructed pCO2 w/ full chemical leverage (totl)
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


%% Calculations -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size;  

st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii dim 

    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    % Set file path
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
        
    % Load and repeat annual value over 12 months of each year (yint)
    co2_am    = reshape(repmat(reshape(ncread([fpath  '/co2_am.nc'],'co2_am'),[],1,dim(3)/12),1,12,1),dim(1),dim(2),[]);

    % Load data
    co2_bp_nl = ncread([fpath  '/co2_bp_nl.nc'],'co2_bp_nl');
    co2_t_nl  = ncread([fpath  '/co2_t_nl.nc'],'co2_t_nl');
    co2_bp    = ncread([fpath  '/co2_bp.nc'],'co2_bp');
    co2_t     = ncread([fpath  '/co2_t.nc'],'co2_t');
        
    % Total pCO2 w/ no chemical levearge
    co2_nl   = co2_bp_nl + co2_t_nl - co2_am;
    
    % Total pCO2 w/ biophysical chemical levearge only
    co2_bpl  = co2_bp + co2_t_nl - co2_am;
    
    % Total pCO2 w/ thermal chemical leverage only
    co2_tl   = co2_bp_nl + co2_t - co2_am;
    
    % Total pCO2 w/ full chemical leverage
    co2_totl = co2_bp + co2_t - co2_am;

    % Save variables to .nc files 
    disp('Saving')     

    nccreate([fpath   '/co2_nl.nc'],'co2_nl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_nl.nc'],'co2_nl',co2_nl)
    ncwriteatt([fpath '/co2_nl.nc'],'co2_nl','units','uatm')
    ncwriteatt([fpath '/co2_nl.nc'],'co2_nl','long_name','recomputed ocean total pCO2 w/o chemical leverage')

    nccreate([fpath   '/co2_bpl.nc'],'co2_bpl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_bpl.nc'],'co2_bpl',co2_bpl)
    ncwriteatt([fpath '/co2_bpl.nc'],'co2_bpl','units','uatm')
    ncwriteatt([fpath '/co2_bpl.nc'],'co2_bpl','long_name','recomputed ocean total pCO2 w/ biophysical chemical leverage')

    nccreate([fpath   '/co2_tl.nc'],'co2_tl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_tl.nc'],'co2_tl',co2_tl)
    ncwriteatt([fpath '/co2_tl.nc'],'co2_tl','units','uatm')
    ncwriteatt([fpath '/co2_tl.nc'],'co2_tl','long_name','recomputed ocean total pCO2 w/ thermal chemical leverage')

    nccreate([fpath   '/co2_totl.nc'],'co2_totl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/co2_totl.nc'],'co2_totl',co2_totl)
    ncwriteatt([fpath '/co2_totl.nc'],'co2_totl','units','uatm')
    ncwriteatt([fpath '/co2_totl.nc'],'co2_totl','long_name','recomputed ocean total pCO2 w/ chemical leverage')
    
    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)  
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~4 min per ensemble member


%% Clean up -----------------------------------

clear;close all;clc
