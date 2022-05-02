%% Figure S1 Calculations 
%
% DESCRIPTION:
%   Carbonate system internal constency calculations.
%   Use ensemble member 1 only
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


%% Load Data and Perform Calculations-----------------------------------

time = ncread([data_path '/Ensemble_Members/time.nc'],'time'); 

fpath = [data_path '/Ensemble_Members/m1/'];
vinfo = ncinfo([fpath 'ta.nc'],'ta');
dim = vinfo.Size;  

% Pre-allocate
co2_calc = NaN(dim(1),dim(2),length(time));

% Perform calculations globally in 5 year chunks
interval = 12*5;
yind     = [1:interval:length(time)]; 

st_0 = now;
for i = 1:length(yind)
    clearvars -except i co2_calc dim st_0 fpath time interval yind 

    disp(['Start 5-year chunk: ' num2str(i) '/' num2str(length(yind))])
    st_1 = now;
    
    if i == length(yind)
        % length time (1812) is not evenly divisible by 5, so change interval for last chunk
        interval = length(time) - yind(i) + 1;
    end   

    % Load 5-year chunk of .nc file and reshape for CO2SYS calculations
    TA   = reshape(ncread([fpath  'ta.nc'],  'ta',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
    DIC  = reshape(ncread([fpath 'dic.nc'],  'dic', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SSS  = reshape(ncread([fpath 'sss.nc'],  'sss', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SST  = reshape(ncread([fpath 'sst.nc'],  'sst', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SIO4 = reshape(ncread([fpath 'sio4.nc'],'sio4', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    PO4  = reshape(ncread([fpath 'po4.nc'],  'po4', [1 1 yind(i)],[Inf Inf interval]), [], 1);
        
    % Pre-allocate blank vector for computed values
    co2Array = NaN(length(TA), 1);
         
    % Find non-NaN values
    xx = find(~isnan(TA) & ~isnan(SSS) & ~isnan(SIO4));
    
    % Perform CO2SYS calculations
    A = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
    
    % Put subsetted values into blank vectors
    co2Array(xx) = A(:, 19);
    
    % Reshape
    co2_calc(:,:,yind(i):yind(i)+interval-1) = reshape(co2Array,dim(1),dim(2),interval);
        
    disp('Finish Calculations')
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)

% Save variable to .nc file 
disp('Saving') 
nccreate([fpath   'co2_calc.nc'],'co2_calc','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
ncwrite( [fpath   'co2_calc.nc'],'co2_calc',co2_calc)
ncwriteatt([fpath 'co2_calc.nc'],'co2_calc','units','uatm')
ncwriteatt([fpath 'co2_calc.nc'],'co2_calc','long_name','pCO2 calculated from model ouptut of DIC and TA')


%% Clean up -----------------------------------

clear;close all;clc
