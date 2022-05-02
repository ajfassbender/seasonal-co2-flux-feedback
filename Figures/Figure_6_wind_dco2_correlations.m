%% Figure 6 Calculations
%
% DESCRIPTION:
%   Local seasonal wind speed vs. dpCO2 correlations.
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
 

%% Make Caluclations -----------------------------------

fpath = [data_path '/Ensemble_Members'];
lat   = ncread([fpath '/lat.nc'],  'lat');   
lon   = ncread([fpath '/lon.nc'],  'lon'); 
time  = ncread([fpath '/time.nc'], 'time'); 

yr   = str2num(datestr(time,'yyyy'));
yri  = find(2090 <= yr & yr < 2100);
mo   = str2num(datestr(time(yri),'mm'));

st_0 = now;
for i = 1:30
    clearvars -except data_path yr yri  mo i st_0 delta delta_nl windy

    disp(['Start ensemble member ' num2str(i)])   
    st_1 = now; 
    
    fpath = [data_path '/Ensemble_Members/m' num2str(i)];

    co2_atm = ncread([fpath  '/co2_atm.nc'],'co2_atm', [1 1 yri(1)], [Inf Inf length(yri)]);
    wind    = ncread([fpath  '/wind.nc'],'wind', [1 1 yri(1)], [Inf Inf length(yri)]);    
    co2_nl  = ncread([fpath  '/co2_corr_nl.nc'],'co2_corr_nl', [1 1 yri(1)], [Inf Inf length(yri)]);
    co2     = ncread([fpath  '/co2_corr.nc'],'co2_corr', [1 1 yri(1)], [Inf Inf length(yri)]);
    
    % 2090s average annual cycles          
    co2_av     = NaN(360,180,12);
    co2_nl_av  = NaN(360,180,12);
    co2_atm_av = NaN(360,180,12);
    wind_av    = NaN(360,180,12);
    for n = 1:12
        clear q;q         = find(mo == n);
        co2_av(:,:,n)     = nanmean(co2(:,:,q),3);
        co2_nl_av(:,:,n)  = nanmean(co2_nl(:,:,q),3);
        co2_atm_av(:,:,n) = nanmean(co2_atm(:,:,q),3);
        wind_av(:,:,n)    = nanmean(wind(:,:,q),3);
    end
    
    delta(:,:,:,i)     = co2_av    - co2_atm_av;
    delta_nl(:,:,:,i)  = co2_nl_av - co2_atm_av;
    windy(:,:,:,i)      = wind_av;
       
    disp(['Complete Ensemble Member ' num2str(i)])
    minutes((now - st_1)*24*60)
end
disp('Finish Calculations')
minutes((now - st_0)*24*60)


%% Wind Speed dpCO2 Correlations -----------------------------------

% Use ensemble mean seasonal cycles

clear wind 
dco2      = nanmean(delta,4);
dco2_nl   = nanmean(delta_nl,4);
wind      = nanmean(windy,4);

R        = NaN.*ones(360,180);
Pval     = NaN.*ones(360,180);
R_nl     = NaN.*ones(360,180);
Pval_nl  = NaN.*ones(360,180);

% Compute correlation for each model grid 
for ii = 1:360
    for j = 1:180
        clear x xtpi y
        x    = squeeze(dco2(ii,j,:));
        x_nl = squeeze(dco2_nl(ii,j,:));
        y    = squeeze(wind(ii,j,:));
        if nansum(x) == 0
        else
            clear r p
            [r,p] = corrcoef(x,y);
            R(ii,j)    = r(2);
            Pval(ii,j) = p(2);
            
            clear r p
            [r,p] = corrcoef(x_nl,y);
            R_nl(ii,j)    = r(2);
            Pval_nl(ii,j) = p(2);
        end
    end
end


%% Save variables to .nc files -----------------------------------
disp('Saving') 

fpath = [data_path '/Ensemble_Members'];
[n,m] = size(R);

vars = {'R','Pval','R_nl','Pval_nl'};
for i = 1:length(vars)
    uvar = char(vars(i));
    eval(['x = ' uvar ';'])
    nccreate([fpath   '/Figure6_wind_dco2_correlations.nc'],uvar,'Dimensions',{'lon',n,'lat',m})
    ncwrite([fpath    '/Figure6_wind_dco2_correlations.nc'],uvar,x)
    ncwriteatt([fpath '/Figure6_wind_dco2_correlations.nc'],uvar,'units','none')
    ncwriteatt([fpath '/Figure6_wind_dco2_correlations.nc'],uvar,'long_name','correlation coefficient or p value')
end


%% Clean Up -----------------------------------

clear;close all;clc
