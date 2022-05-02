%% step0_nc_reprocess.m
%
% DESCRIPTION:
% Reprocess model files: update variable names and units. 
%   Save reprocessed data into new .nc files.
%   This script assumes that ensemble member data files are grouped by variable and
%       saved in subfolders titled with the variable names. These folders reside within an 
%       'Ensemble_Members' subfolder in the 'data_path' directory: 
%       [data_path '/Ensemle_Members/*VariableName*/*files*']
%   This script assumes that control run files are located in a subfolder titled
%        'Control_Run' in the 'data_path' directory: [data_path '/Control_Run/*files*']
%   Required ESM2M surface fields, previously regridded to 1° latitide by 1° longitude, can be
%   accessed here: http://poseidon.princeton.edu/
%
%
% USER INPUT:
%   data_path: directory containing 'Ensemble_Members' and 'Control_Run' subfolders
%   save_path: directory where reprocessed .nc files will be saved
%
% OUTPUT: 
%   ensemble member .nc files: 'ta','dic','sss','sst','dco2','mflux','ice_frac','mld','co2',
%                              'co2_atm','po4','sio4','wind'
%   control run .nc files:     'ta','dic','sss','sst','dco2','mflux','mld','co2', 'po4','sio4'
%
% FUNCTIONS CALLED:
%   daynoleap2datenum: https://github.com/kakearney/daynoleap2datenum-pkg
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

% Directory containing 'Ensemble_Members' and 'Control_Run' subfolders
disp('set your path')
data_path = cd;

% Directory where reprocessed .nc files will be saved 
disp('set your path')
save_path = cd;


%% Reprocess Transient Ensemble Member Files -----------------------------------
             
% Existing folder and variable names and the new variable names for reprocessed .nc files
var_model(:,1)   = {'SFC_ALK','SFC_DIC','SSS','SST','dic_deltap','dic_kw','dic_stf','ice_frac','mld','pCO2','PCO2_ATM','sfc_po4','sfc_sio4','wind'};
var_rename(:,1)  = {'ta','dic','sss','sst','dco2','k_pv','mflux','ice_frac','mld','co2','co2_atm', 'po4','sio4','wind'};
model_units(:,1) = {'mol/kg','mol/kg','none','K','uatm','m/s','mole/m^2/sec','none','m','uatm','atm','mol/kg','mol/kg','m/s'};
final_units(:,1) = {'umol/kg','umol/kg','none','C','uatm','m/s','mol/m2/d','none','m','uatm','atm','umol/kg','umol/kg','m/s'};
scale_factor     = [10^6;10^6;1;-273.15;1;1;86400;1;1;1;1;10^6;10^6;1];
t1 = cell2table([var_model model_units var_rename final_units],'VariableNames',{'Model Var','Model Units','New Var','New Units'});
t2 = [t1 table(scale_factor,'VariableNames',{'Scale Factor'})]

% Cycle through variable folders
st_0 = now;
for j = 1:length(var_model)
    st_1 = now;
    var_name     = char(var_model(j));
    var_new_name = char(var_rename(j));
    disp(['Start reprocessing ' var_name])

    % Get list of .nc files in folder
    fpath = [data_path '/Ensemble_Members/' var_name];
    flist = dir(fullfile(fpath,'*.nc'));
    
    % Cylce through ensemble member files
    for i = 1:length(flist)
        st_2 = now;

        % Directory to save reprocessed .nc file
        spath = [save_path '/Ensemble_Members/m' num2str(i)];
        sname = [spath '/' var_new_name '.nc'];

        % Determine if save directory exists
        if ~exist(spath, 'dir')
            % Make directory if needed
            mkdir(spath)
        end

        % Has the variable already been reprocessed for this ensemble member?
        if exist(sname) == 2
        else
            fname = [fpath '/' flist(i).name];

            % Load data
            var = single(ncread(fname,var_name)); % variable
    
            if j == 11 % no units in pCO2_ATMm file
                units = 'uatm';
            else
                units = ncreadatt(fname,var_name,'units');
            end

            info  = ncreadatt(fname,var_name,'long_name');
            vinfo = ncinfo(fname,var_name);
            dim = vinfo.Size;
            fill_val = ncreadatt(fname,var_name,'_FillValue');
    
            % Set missing values to NaN
            z = find(var == fill_val);
            var(z) = NaN;
    
            % Apply scaling if needed
            if scale_factor(j) == 1
            elseif scale_factor(j) == -273.15 
                var = var + scale_factor(j); % addition: K to C
            else
                var = var.* scale_factor(j); % multiplication: unit changes
            end
    
            % SSS and SST files are the wrong size
            if j > 2 && j < 5
                var = squeeze(var);
                dim = size(var);
            end
            
            % Save variable w/ attribute information to .nc file
            nccreate(sname,var_new_name,'Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
            ncwrite(sname,var_new_name,var)
            ncwriteatt(sname,var_new_name,'units',char(final_units(j)))
            ncwriteatt(sname,var_new_name,'long_name',info)
        end
        disp(['Finish ensemble member ' num2str(i)])
        minutes((now - st_2)*24*60)   
    end
    disp(['Finish ' var_name])
    minutes((now - st_1)*24*60)
end
disp(['Complete reprocessing'])
minutes((now - st_0)*24*60)

%% Save time, lat, and lon in Ensemble_Members folder -----------------------------------

fpath = [data_path '/Ensemble_Members/SSS'];
flist = dir(fullfile(fpath,'*.nc'));
fname = [fpath '/' flist(1).name];

% Load lon, lat, and time
lon    = single(ncread(fname,'xt_ocean'));  
lat    = single(ncread(fname,'yt_ocean'));
ttime  = single(ncread(fname,'time'));

% Save lon and lat to .nc file
spath = [save_path '/Ensemble_Members'];
dim = size(lon);
nccreate([spath '/lon.nc'],'lon','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/lon.nc'],'lon',lon)
dim = size(lat);
nccreate([spath '/lat.nc'],'lat','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/lat.nc'],'lat',lat)

% Save time to .nc file (correct for no leap years in the model)
time = daynoleap2datenum(ttime, 1950, 'dn');
dim = size(time);
nccreate([spath '/time.nc'],'time','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/time.nc'],'time',time)


%% Reprocess Control Run Files -----------------------------------

clearvars -except data_path save_path

% Existing variable (and folder) names and the new names for reprocessed .nc files
var_model(:,1)   = {'sfc_alk','sfc_dic','sss','sst','dic_deltap','dic_stf','mld','pco2surf','sfc_po4','sfc_sio4'};
var_rename(:,1)  = {'ta','dic','sss','sst','dco2','mflux','mld','co2', 'po4','sio4'};
model_units(:,1) = {'mol/kg','mol/kg','none','K','uatm','mole/m^2/sec','m','uatm','mol/kg','mol/kg'};
final_units(:,1) = {'umol/kg','umol/kg','none','C','uatm','mol/m2/d','m','uatm','umol/kg','umol/kg'};
scale_factor     = [10^6;10^6;1;-273.15;1;86400;1;1;10^6;10^6];
t1 = cell2table([var_model model_units var_rename final_units],'VariableNames',{'Model Var','Model Units','New Var','New Units'});
t2 = [t1 table(scale_factor,'VariableNames',{'Scale Factor'})]

% Get list of variable file names from control_run folder
flist = dir(fullfile([data_path '/Control_Run'],'*.nc'));

% Determine if save directory exists
spath = [save_path '/Control_Run'];
if ~exist(spath, 'dir')
    % Make directory if needed
    mkdir(spath)
end

% Cycle through variables and reprocess
st_0 = now;
for i = 1:length(flist)
    st_1 = now;

    % Get variable name from file name (format: variable_1x1_1861_2100.nc)
    file_name    = char(flist(i).name);
    vinds        = strfind(file_name, '_1x1'); 
    var_name     = char(file_name(1:vinds(1)-1));
    var_ind      = strcmp(var_name,var_model);
    var_new_name = char(var_rename(var_ind));
    sname        = [spath '/' var_new_name '.nc'];

    % Has the variable already been reprocessed for this ensemble member?
    if exist(sname) == 2 || nansum(var_ind) == 0
    else
        disp(['Start reprocessing ' var_name])

        % Load variable
        full_path = [data_path '/Control_Run/' file_name]; 
        var   = ncread(full_path,var_name); 
        units = ncreadatt(full_path,var_name,'units');
        info  = ncreadatt(full_path,var_name,'long_name');
        vinfo = ncinfo(full_path,var_name);
        dim   = vinfo.Size;
        missing_val = ncreadatt(full_path,var_name,'missing_value');
    
        % Set missing values to NaN
        z = find(var == missing_val);
        var(z) = NaN;
    
        % Apply scaling if needed
        if scale_factor(var_ind) == 1
        elseif scale_factor(var_ind) == -273.15 
            var = var + scale_factor(var_ind); % addition: K to C
        else
            var = var.* scale_factor(var_ind); % multiplication: unit changes
        end
    
        % Save variable w/ attribute information to .nc file
        nccreate(sname,var_new_name,'Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
        ncwrite(sname,var_new_name,var)
        ncwriteatt(sname,var_new_name,'units',char(final_units(var_ind)))
        ncwriteatt(sname,var_new_name,'long_name',info)
    end
    disp(['Complete'])
    minutes((now - st_1)*24*60)
end
disp(['Total time'])
minutes((now - st_0)*24*60)
% ~10 min total

%% Save time, lat, and lon in Control_Run folder -----------------------------------

fname = [data_path '/Control_Run/dic_deltap_1x1_1861_2100.nc'];

% Load lon, lat, and time
lon    = single(ncread(fname,'xt_ocean'));  
lat    = single(ncread(fname,'yt_ocean'));
ttime  = single(ncread(fname,'time'));

% Save lon and lat to .nc file
spath = [save_path '/Control_Run'];
dim = size(lon);
nccreate([spath '/lon.nc'],'lon','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/lon.nc'],'lon',lon)
dim = size(lat);
nccreate([spath '/lat.nc'],'lat','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/lat.nc'],'lat',lat)

% Save time to .nc file (correct for no leap years in the model)
time = daynoleap2datenum(ttime, 0001, 'dn');
dim = size(time);
nccreate([spath '/time.nc'],'time','Dimensions',{'r',dim(1),'c',dim(2)},'DeflateLevel',5)
ncwrite([spath '/time.nc'],'time',time)


%% Clean up

clear;close all;clc
