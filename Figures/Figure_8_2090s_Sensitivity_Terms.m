%% Figure 8 Calculations 
%
% DESCRIPTION:
%   Calculate 2090s ensemble mean and control run sensitivity terms.
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

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size;  

% Pre-Allocate
rf_2090       = NaN(dim(1),dim(2),30);
af_2090       = NaN(dim(1),dim(2),30);
co2_am_2090   = NaN(dim(1),dim(2),30);

st_0 = now;
for i = 1:30
    clearvars -except data_path i yri uyri st_0 rf_2090 af_2090 co2_am_2090

    disp(['Start ensemble member ' num2str(i)])   
    st_1 = now; 
    
    fpath = [data_path '/Ensemble_Members/m' num2str(i)];

     rf_2090(:,:,i)      = nanmean(ncread([fpath  '/rf.nc'],'rf',[1 1 yri(1)], [Inf Inf length(yri)]),3);
     af_2090(:,:,i)      = nanmean(ncread([fpath  '/af.nc'],'af',[1 1 yri(1)], [Inf Inf length(yri)]),3);
     co2_am_2090(:,:,i)  = nanmean(ncread([fpath  '/co2_am.nc'],'co2_am',[1 1 uyri(1)], [Inf Inf length(uyri)]),3);

   disp(['Complete Ensemble Member ' num2str(i)])
   minutes((now - st_1)*24*60)
end
disp('Finish Calculations')
minutes((now - st_0)*24*60)


%% Compute Ensemble Mean Values -----------------------------------

rf      = nanmean(rf_2090,3);
af      = nanmean(af_2090,3);
co2_am  = nanmean(co2_am_2090,3);


%% Compute Control Run Values
disp('Control Run Calculations')
fpath = [data_path '/Control_Run'];

rf_pi      = nanmean(ncread([fpath  '/rf.nc'],'rf'),3);
af_pi      = nanmean(ncread([fpath  '/climatology.nc'],'af_mean'),3);
co2_am_pi  = nanmean(ncread([fpath  '/co2_am.nc'],'co2_am'),3);

%% Save variables to .nc files -----------------------------------
disp('Saving')

fpath = [data_path '/Ensemble_Members'];
[n,m] = size(co2_am);

vars  = {'rf','af','co2_am','rf_pi','af_pi','co2_am_pi'};
uni   = {'none','none','uatm','none','none','uatm'};
for i = 1:length(vars)
    uvar = char(vars(i));
    eval(['x = ' uvar ';'])
    nccreate([fpath   '/Figure8_2090s_Sensitivity_Terms.nc'],uvar,'Dimensions',{'lon',n,'lat',m})
    ncwrite([fpath    '/Figure8_2090s_Sensitivity_Terms.nc'],uvar,x)
    ncwriteatt([fpath '/Figure8_2090s_Sensitivity_Terms.nc'],uvar,'units',char(uni(i)))
    ncwriteatt([fpath '/Figure8_2090s_Sensitivity_Terms.nc'],uvar,'long_name','2090s or Control Run (PI) annual mean or climatological values')
end


%% -------- Clean Up

clear;close all;clc
