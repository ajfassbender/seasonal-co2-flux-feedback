%% Figure 4 Calculations
%
% DESCRIPTION:
%   1950s to 2090s ensemble mean pCO2 amplitude changes.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;
fig_path  = cd;

%% Calculate 1950s and 2090s Amplitudes -----------------------------------

lat  = ncread([data_path '/Ensemble_Members/lat.nc'],  'lat');   
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],  'lon');  
time = ncread([data_path '/Ensemble_Members/time.nc'],'time'); 

yr      = str2num(datestr(time,'yyyy'));
uyr     = unique(yr);
ind_50s = find(yr < 1960);
ind_90s = find(2090 <= yr & yr < 2100);

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim   = vinfo.Size; 

vars = {'dic','ta','sst','wind','co2_totl','mld','co2_bp','co2_t','co2_t_nl','co2_bp_nl'};
for i = 1:length(vars)
    uvar = char(vars(i));
    if strcmp(uvar,'mld') == 1        
        eval([uvar '_max_50s = NaN(dim(1),dim(2),30);']);
        eval([uvar '_max_90s = NaN(dim(1),dim(2),30);']);
    else
        eval(['A' uvar '_50s = NaN(dim(1),dim(2),30);']);
        eval(['A' uvar '_90s = NaN(dim(1),dim(2),30);']);
    end
end

st_0 = now;
for ii = 1:30
    st_1 = now;

    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    disp(['Start ensemble member ' num2str(ii)])

    % Calculate values for the two decades
    for j = 1:2
        disp(['Start processing decade ' num2str(j)])   
        if j == 1
            ind = ind_50s;
            txt = '_50s';
        else
            ind = ind_90s;
            txt = '_90s';
        end 
    
        for i = 1:length(vars)
            uvar = char(vars(i));
            x    = ncread([fpath '/' uvar '.nc'],uvar,[1 1 ind(1)],[Inf Inf length(ind)]);
            if strcmp(uvar,'mld') == 1
                eval([uvar '_max' txt '(:,:,ii) = nanmean(squeeze(max(reshape(x,dim(1),dim(2),12,[]),[],3)),3);']);
            else
                x = eval(uvar);
                eval(['A' uvar txt '(:,:,ii) = nanmean(squeeze(max(reshape(x,dim(1),dim(2),12,[]),[],3) - min(reshape(x,dim(1),dim(2),12,[]),[],3)),3);']);
            end
            clear x uvar
        end
    end
    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Chemical Leverage Impact on pCO2 BP and pCO2 T

ilat       = -89.5:1:89.5;
ilon       =  0.5:1:359.5;

dco2_bp    = ncread(fpath,'Aco2_bp_90s')    - ncread(fpath,'Aco2_bp_50s');
dco2_bp_nl = ncread(fpath,'Aco2_bp_nl_90s') - ncread(fpath,'Aco2_bp_nl_50s');
dco2_bp_p  = 100 .* ((dco2_bp-dco2_bp_nl)./dco2_bp);
area_weighted_mean(dco2_bp_p,ilat,ilon) % 45 ± 88 percent
area_weighted_std(dco2_bp_p,ilat,ilon)

dco2_t    = ncread(fpath,'Aco2_t_90s')    - ncread(fpath,'Aco2_t_50s');
dco2_t_nl = ncread(fpath,'Aco2_t_nl_90s') - ncread(fpath,'Aco2_t_nl_50s');
dco2_t_p  = 100 .* ((dco2_t-dco2_t_nl)./dco2_t);
area_weighted_mean(dco2_t_p,ilat,ilon) % 101 ± 11 percent
area_weighted_std(dco2_t_p,ilat,ilon) 


% Save Results -----------------------------------

fpath = [data_path '/Ensemble_Members'];

vars = {'dic','ta','sst','wind','co2_totl','mld','co2_bp','co2_t','co2_t_nl','co2_bp_nl'};
uni  = {'umol/kg','umol/kg','C','m/s','uatm','m','uatm','uatm','uatm','uatm'};

for i = 1:length(vars)
    if strcmp(vars(i),'mld') == 1
        uvar = [char(vars(i)) '_max_50s'];
        eval(['x = ' uvar ';']);
        X = nanmean(x,3);
        nccreate([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'Dimensions',{'lon',dim(1),'lat',dim(2)})
        ncwrite( [fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,X)
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'units',char(uni(i)))
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'long_name','1950s mean max annual MLD')

        uvar = [char(vars(i)) '_max_90s'];
        eval(['x = ' uvar ';']);
        X = nanmean(x,3);
        nccreate([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'Dimensions',{'lon',dim(1),'lat',dim(2)})
        ncwrite( [fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,X)
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'units',char(uni(i)))
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'long_name','2090s mean max annual MLD')
    else
        uvar = ['A' char(vars(i)) '_50s'];
        eval(['x = ' uvar ';']);
        X = nanmean(x,3);
        nccreate([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'Dimensions',{'lon',dim(1),'lat',dim(2)})
        ncwrite( [fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,X)
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'units',char(uni(i)))
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'long_name','1950s mean seasonal amplitude')

        uvar = ['A' char(vars(i)) '_90s'];
        eval(['x = ' uvar ';']);
        X = nanmean(x,3);
        nccreate([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'Dimensions',{'lon',dim(1),'lat',dim(2)})
        ncwrite( [fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,X)
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'units',char(uni(i)))
        ncwriteatt([fpath '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'],uvar,'long_name','2090s mean sesaaonal amplitude')
    end
end 


%% Area-weighted changes globally -----------------------------------

clearvars -except data_path fig_path

vars = {'Adic','Ata','Asst','Awind','mld_max','Aco2_totl','Aco2_t','Aco2_t_nl','Aco2_bp','Aco2_bp_nl'};
uni  = {'umol/kg','umol/kg','C','m/s','m','uatm','uatm','uatm','uatm','uatm'};

fpath = [data_path '/Ensemble_Members'];
area  = reshape(transpose(getArea(ilat,ilon)),[],1);

for i = 1:length(vars)
    X_50 = ncread([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'], [char(vars(i)) '_50s']);  
    X_90 = ncread([fpath   '/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'], [char(vars(i)) '_90s']);   
    
    delta = X_90-X_50;
    wmean(i,1)  = area_weighted_mean(delta,ilat,ilon);
    wstdev(i,1) = area_weighted_std( delta,ilat,ilon);

    delta_p = 100 .* delta./X_50;
    wmean_p(i,1)  = area_weighted_mean(delta_p,ilat,ilon);
    wstdev_p(i,1) = area_weighted_std( delta_p,ilat,ilon);
end 

t0 = cell2table(vars(:),'VariableNames',{'1950s-2090s Change'});
t1 = cell2table(uni(:),'VariableNames',{'Units'});
t2 = array2table([wmean wstdev wmean_p wstdev_p],'VariableNames',{'Weighted Mean' 'Weighted StDev' 'Weighted Mean (%)' 'Weighted StDev (%)'});
t  = [t0 t2 t1]
save([fig_path '/Ensemble_Mean_AmplitudeChanges_1950_2090.mat'], 't')

%% ----- Clean up

%clear;close all;clc
