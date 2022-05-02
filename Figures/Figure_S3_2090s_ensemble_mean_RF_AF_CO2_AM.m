%% Figure S3 Calculations 
%
% DESCRIPTION:
%   pCO2 AM, RF, and AF changes from the pre-industrial to the 2090s.
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
fig_path  = cd;


%% Ensemble Mean 2090s Values and StDevs -----------------------------------

fpath = [data_path '/Ensemble_Members/'];
vinfo = ncinfo([fpath 'm1/ta.nc'],'ta');
dim   = vinfo.Size;  

time = ncread([fpath 'time.nc'],'time'); 
yr   = str2num(datestr(time,'yyyy'));
uyr  = unique(yr);
yri  = find(2090 <= yr  &  yr < 2100);
uyri = find(2090 <= uyr & uyr < 2100);

% Pre-allocate
RF = NaN(dim(1),dim(2),30);
AF = NaN(dim(1),dim(2),30);
AM = NaN(dim(1),dim(2),30);
RF_sd = NaN(dim(1),dim(2),30);
AF_sd = NaN(dim(1),dim(2),30);
AM_sd = NaN(dim(1),dim(2),30);

st_0 = now;
for ii = 1:30
    clear co2_am rf af

    st_1 = now;

    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    disp(['Start ensemble member ' num2str(ii)])
    
    % Load Data
    rf      = ncread([fpath '/rf.nc'],'rf',[1 1 yri(1)],[Inf Inf length(yri)]);
    af      = ncread([fpath '/af.nc'],'af',[1 1 yri(1)],[Inf Inf length(yri)]);
    co2_am  = ncread([fpath '/co2_am.nc'],'co2_am',[1 1 uyri(1)],[Inf Inf length(uyri)]);
    
    RF(:,:,ii)     = nanmean(rf,3);
    AF(:,:,ii)     = nanmean(af,3);
    AM(:,:,ii)     = nanmean(co2_am,3);
    RF_sd(:,:,ii)  = nanstd(rf,[],3);
    AF_sd(:,:,ii)  = nanstd(af,[],3);
    AM_sd(:,:,ii)  = nanstd(co2_am,[],3);

    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)

rf_90s = nanmean(RF,3);
af_90s = nanmean(AF,3);
co2_am_90s = nanmean(AM,3);
rf_90s_sd = nanmean(RF_sd,3);
af_90s_sd = nanmean(AF_sd,3);
co2_am_90s_sd = nanmean(AM_sd,3);


%% Save variables to .nc files -----------------------------------
disp('Saving')

fpath = [data_path '/Ensemble_Members'];
[n,m] = size(rf_90s);

vars  = {'rf_90s','af_90s','co2_am_90s','rf_90s_sd','af_90s_sd','co2_am_90s_sd'};
uni   = {'none','none','uatm','none','none','uatm'};
for i = 1:length(vars)
    uvar = char(vars(i));
    eval(['x = ' uvar ';'])
    nccreate([fpath   '/FigureS3_2090s_Ensemble_Mean_RF_AF_CO2_AM.nc'],uvar,'Dimensions',{'lon',n,'lat',m})
    ncwrite([fpath    '/FigureS3_2090s_Ensemble_Mean_RF_AF_CO2_AM.nc'],uvar,x)
    ncwriteatt([fpath '/FigureS3_2090s_Ensemble_Mean_RF_AF_CO2_AM.nc'],uvar,'units',char(uni(i)))
    ncwriteatt([fpath '/FigureS3_2090s_Ensemble_Mean_RF_AF_CO2_AM.nc'],uvar,'long_name','2090s ensemlbe mean values')
end

%% Clean up -----------------------------------

clear;close all;clc
