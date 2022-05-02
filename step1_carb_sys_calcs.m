%% step1_carb_sys_calcs.m
%
% DESCRIPTION:
%    Calculate DIC and Revelle Factor values for each ensemble member and the control run.
%
% INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT: 
%   rf         (unitless) - for control run and ensemble members
%   dic_calc   (micromole per kg) - for control run and ensemble members
%   climatolgy (multiple) - file containing 12-month climatologies and standard deviations
%                           for control run variables: 'co2','dic','ta','dic_calc','rf','dco2',
%                           'mflux','mld','po4','sio4','sss','sst'
%
% FUNCTIONS CALLED:
%   CO2SYS: Using dissociation constants of Mehrbach et al. (1973), as refitted by Dickson 
%   and Millero (1987), the hydrogen sulfate dissociation constant of Dickson et al. (1990),
%   and the boron-to-chlorinity ratio of Uppström 1974.
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your path')
data_path = cd;

%% Control Run Carbonate System Calculations -----------------------------------

fpath = [data_path '/Control_Run/'];
time  = ncread([fpath '/time.nc'],'time');

% Get final file size from existing variable (ta)
vinfo = ncinfo([fpath 'ta.nc'],'ta');
dim   = vinfo.Size; 

% Perform calculations globally in 5 year chunks
interval = 12*5;
yind = [1:interval:length(time)];  
st_0 = now;   
for i = 1 : length(yind)
    disp(['Start 5-year chunk: ' num2str(i) '/' num2str(length(yind))])
    st_1 = now;

    clear rfArray dicArray xx A rf dic_calc TA CO2 SSS SST SIO4 PO4 dim
    
    % Load 5-year chunk of .nc file and reshape for CO2SYS calculations
    TA   = reshape(ncread([fpath   'ta.nc'],'ta',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
    CO2  = reshape(ncread([fpath  'co2.nc'],'co2', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SSS  = reshape(ncread([fpath  'sss.nc'],'sss', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SST  = reshape(ncread([fpath  'sst.nc'],'sst', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    SIO4 = reshape(ncread([fpath 'sio4.nc'],'sio4',[1 1 yind(i)],[Inf Inf interval]), [], 1);
    PO4  = reshape(ncread([fpath  'po4.nc'],'po4', [1 1 yind(i)],[Inf Inf interval]), [], 1);
    
    % Pre-allocate blank vectors for computed values
    rfArray   = NaN(length(TA), 1);
    dicArray  = NaN(length(TA), 1);
     
    % Find non-NaN values
    xx = find(~isnan(TA) & ~isnan(SSS) & ~isnan(SIO4));

    % Perform CO2SYS calculations
    A = CO2SYS(TA(xx),CO2(xx),1,4,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
    
    % Add computed values to blank vectors
    rfArray(xx)   = A(:, 29);
    dicArray(xx)  = A(:,  2);

    % Reshape
    rf       = reshape(rfArray,  dim(1), dim(2), interval);
    dic_calc = reshape(dicArray, dim(1), dim(2), interval);
    
    % Save variables to .nc files 
    disp('Saving')
    if i == 1;
        % Create new .nc file of correct dimensions          
        nccreate([fpath 'rf.nc'],'rf','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
        
        % Write 5-year variable chunk and attributes
        ncwrite([fpath 'rf.nc'],'rf',rf,[1 1 yind(i)])
        ncwriteatt([fpath 'rf.nc'],'rf','units','none')
        ncwriteatt([fpath 'rf.nc'],'rf','long_name','Revelle Factor')

        % Create new .nc file of correct dimensions  
        nccreate([fpath 'dic_calc.nc'],'dic_calc','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
       
        % Write 5-year variable chunk and attributes
        ncwrite([fpath 'dic_calc.nc'],'dic_calc',dic_calc,[1 1 yind(i)])
        ncwriteatt([fpath 'dic_calc.nc'],'dic_calc','units','umol/kg')
        ncwriteatt([fpath 'dic_calc.nc'],'dic_calc','long_name','DIC calculated from model TA and pCO2')
    else
        % Write 5-year variable chunk
        ncwrite([fpath 'rf.nc'],'rf',rf,[1 1 yind(i)])
        ncwrite([fpath 'dic_calc.nc'],'dic_calc',dic_calc,[1 1 yind(i)])
    end
    disp(['Finish 5-year chunk: ' num2str(i) '/' num2str(length(yind))])
    minutes((now - st_1)*24*60)
end
disp('Complete Control Run CO2Sys Calcs')
minutes((now - st_0)*24*60)
 ~40 minutes

%% Control Run 12-Month Climatologies -----------------------------------

clearvars -except base_path data_path
clc

% Get file names for variables to process
fpath = [data_path '/Control_Run/'];
flist = dir(fullfile([data_path '/Control_Run'],'*.nc'));
file_name(:,1) = {flist.name};

% Get index of files names excluding time, lon, or lat
x = find(strcmp(file_name,'lon.nc')==0 & strcmp(file_name,'lat.nc')==0 & strcmp(file_name,'time.nc')==0);
fname = file_name(x);

% Cycle through variables to make 12-month climatologies
st_0 = now;
for i = 1:length(fname)
    st_1 = now;
    clear var_mean var_stdev var_mn var_st

    % Load Variable
    fvar      = char(fname(i));
    v_name    = fvar(1:end-3);
    var       = ncread([fpath fvar],v_name);
    var_units = ncreadatt([fpath fvar],v_name,'units');

    % Reshape variable
    [n,m,o] = size(var);
    VAR     = reshape(var,[],o);

    for j = 1:12      
        % Compute average of all monthly values  
        var_mn = nanmean(VAR(:,j:12:end),2);
        var_mean(:,:,j) = reshape(var_mn,n,m);

        % Compute standard deviation of all monthly values
        var_st = nanstd(VAR(:,j:12:end),[],2);
        var_stdev(:,:,j) = reshape(var_st,n,m);
    end

    % Save 
    disp('Saving')
    nccreate([fpath 'climatology.nc'],[v_name '_mean'],'Dimensions',{'lon',n,'lat',m,'time' 12},'DeflateLevel',5)
    ncwrite( [fpath 'climatology.nc'],[v_name '_mean'],var_mean)
    ncwriteatt([fpath 'climatology.nc'],[v_name '_mean'],'units',var_units)
    ncwriteatt([fpath 'climatology.nc'],[v_name '_mean'],'long_name','Monthly climatology over all years')

    nccreate([fpath 'climatology.nc'],[v_name '_stdev'],'Dimensions',{'lon',n,'lat',m,'time' 12},'DeflateLevel',5)
    ncwrite( [fpath 'climatology.nc'],[v_name '_stdev'],var_stdev)
    ncwriteatt([fpath 'climatology.nc'],[v_name '_stdev'],'units',var_units)
    ncwriteatt([fpath 'climatology.nc'],[v_name '_stdev'],'long_name','Monthly standard deviation over all years')

    disp(['Complete variable ' v_name])
    minutes((now - st_1)*24*60)
end
disp('Complete climatology calcs')
minutes((now - st_0)*24*60)


%% Ensemble Member Carbonate System Calculations -----------------------------------

clearvars -except base_path data_path
clc

fpath = [data_path '/Ensemble_Members/'];
time  = ncread([fpath '/time.nc'],'time');

% Get final file size from existing variable (ta)
vinfo = ncinfo([fpath 'm1/ta.nc'],'ta');
dim = vinfo.Size;  

% Cycle through 30 ensemble members
st_0 = now;
for ii = 1:30
    clearvars -except time ii st_0 fpath data_path dim

    st_1 = now;
    disp(['Start ensemble member ' num2str(ii)])

    interval = 12*5;
    yind = [1:interval:length(time)]; 

    % Pre-allocate
    rf       = NaN(dim);
    dic_calc = NaN(dim);

    % Perform calculations globally in 5 year chunks
    for i = 1 : length(yind)
        disp(['Start 5-year chunk: ' num2str(i) '/' num2str(length(yind))])
        st_1 = now;

        if i == length(yind)
        	% length time (1812) is not evenly divisible by 5, so change interval for last chunk
        	interval = length(time) - yind(i) + 1;
        end
    
        clear rfArray dicArray xx A TA CO2 SSS SST SIO4 PO4
        
        % Load 5-year chunk of .nc file and reshape for CO2SYS calculations
        TA   = reshape(ncread([fpath 'm' num2str(ii)   '/ta.nc'],'ta',   [1 1 yind(i)],[Inf Inf interval]), [], 1);
        CO2  = reshape(ncread([fpath 'm' num2str(ii)  '/co2.nc'],'co2',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
        SSS  = reshape(ncread([fpath 'm' num2str(ii)  '/sss.nc'],'sss',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
        SST  = reshape(ncread([fpath 'm' num2str(ii)  '/sst.nc'],'sst',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
        SIO4 = reshape(ncread([fpath 'm' num2str(ii) '/sio4.nc'],'sio4', [1 1 yind(i)],[Inf Inf interval]), [], 1);
        PO4  = reshape(ncread([fpath 'm' num2str(ii)  '/po4.nc'],'po4',  [1 1 yind(i)],[Inf Inf interval]), [], 1);
        
        % Pre-allocate blank vectors for computed values
        rfArray   = NaN(length(TA), 1);
        dicArray  = NaN(length(TA), 1);
         
        % Find non-NaN values
        xx = find(~isnan(TA) & ~isnan(SSS) & ~isnan(SIO4));
    
        % Perform CO2SYS calculations
        A = CO2SYS(TA(xx),CO2(xx),1,4,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
        
        % Add computed values to blank vectors
        rfArray(xx)   = A(:, 29);
        dicArray(xx)  = A(:,  2);
    
        % Reshape
        rf(:,:,yind(i):yind(i)+interval-1)        = reshape(rfArray,  dim(1), dim(2), interval);
        dic_calc(:,:,yind(i):yind(i)+interval-1)  = reshape(dicArray, dim(1), dim(2), interval);
    end
    disp('Finish Calculations')
    minutes((now - st_1)*24*60)

    % Save variables to .nc files 
    disp('Saving')
       
    nccreate([fpath   'm' num2str(ii) '/rf.nc'],'rf','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   'm' num2str(ii) '/rf.nc'],'rf',rf);
    ncwriteatt([fpath 'm' num2str(ii) '/rf.nc'],'rf','units','none');
    ncwriteatt([fpath 'm' num2str(ii) '/rf.nc'],'rf','long_name','Revelle Factor');

    nccreate([fpath   'm' num2str(ii) '/dic_calc.nc'],'dic_calc','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   'm' num2str(ii) '/dic_calc.nc'],'dic_calc',dic_calc);
    ncwriteatt([fpath 'm' num2str(ii) '/dic_calc.nc'],'dic_calc','units','umol/kg');
    ncwriteatt([fpath 'm' num2str(ii) '/dic_calc.nc'],'dic_calc','long_name','DIC calculated from model TA and pCO2');

    disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~13 min/ensemble member


%% Clean Up -----------------------------------

clear;close all;clc
