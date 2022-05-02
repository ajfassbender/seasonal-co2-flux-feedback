%% step4_pco2bp_calcs.m
%
% DESCRIPTION:
%   Calculate biophysical pCO2 component (pCO2 BP) w/ and w/o chemical leverage.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT:    
%   climatolgy (unitless) - file containing 12-month climatologies and standard deviations
%                           for control run variables: 'sf' and 'af'
%   sf        (unitless) - ensemble member salnity factor 
%   af        (unitless) - ensemble member alkalinity factor
%   co2_bp_nl (microatmosphere) - ensemble member biophysical pCO2 component w/o chemical leverage (nl)
%   co2_bp    (microatmosphere) - ensemble member biophysical pCO2 component w/ chemical leverage
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

%% Control Run 12-Month Alaklinity & Salinity Factor Climatologies -----------------------------------
 
% Set file path
fpath = [data_path '/Control_Run/'];

vinfo = ncinfo([fpath 'climatology.nc'],'ta_mean');
dim   = vinfo.Size;  

% Pre-allocate
af    = NaN(dim(1),dim(2),dim(3));
sf    = NaN(dim(1),dim(2),dim(3));
 
st_0 = now;
for i = 1:12
    clearvars -except data_path st_0 i fpath dim af sf af_stdev sf_stdev
     
    % Load monthly climatology maps and reshape to vectors
    DIC  = reshape(ncread([fpath 'climatology.nc'],'dic_calc_mean',[1 1 i],[Inf Inf 1]), [], 1);
    TA   = reshape(ncread([fpath 'climatology.nc'],'ta_mean',  [1 1 i],[Inf Inf 1]), [], 1);
    SSS  = reshape(ncread([fpath 'climatology.nc'],'sss_mean', [1 1 i],[Inf Inf 1]), [], 1);
    SST  = reshape(ncread([fpath 'climatology.nc'],'sst_mean', [1 1 i],[Inf Inf 1]), [], 1);
    SIO4 = reshape(ncread([fpath 'climatology.nc'],'sio4_mean',[1 1 i],[Inf Inf 1]), [], 1);
    PO4  = reshape(ncread([fpath 'climatology.nc'],'po4_mean', [1 1 i],[Inf Inf 1]), [], 1);
         
    % Find non-NaN values
    xx = find(~isnan(TA) & ~isnan(DIC) & ~isnan(SIO4));
 
    % Perform CO2SYS calculations 
    A = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
    CO2     = NaN(length(TA),1);
    CO2(xx) = A(:,19);
     
    % Define purturbation
    ta_add  = [-0.01 0.01];
    sal_add = [-0.0001 0.0001];
 
    % Pre-allocate
    build_af = NaN(length(TA), length(ta_add));
    build_sf = NaN(length(TA), length(ta_add));

    % Calculate sensitivity terms
    for j = 1:length(ta_add)
        % Alkalinity Factor
        A    = CO2SYS(TA(xx) + ta_add(j),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
        dco2 = A(:,19) - CO2(xx);
        build_af(xx,j) = (dco2 ./ CO2(xx)) .* (ta_add(j) ./ TA(xx)).^-1;   
        clear A dco2
        
        % Salinity Factor
        A    = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx) + sal_add(j) ,SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
        dco2 = A(:,19) - CO2(xx);
        build_sf(xx,j) = (dco2 ./ CO2(xx)) .* (sal_add(j) ./ SSS(xx)).^-1; 
        clear A dco2
    end
 
    af(:,:,i)  = reshape(nanmean(build_af,2),  dim(1),dim(2));
    sf(:,:,i)  = reshape(nanmean(build_sf,2),  dim(1),dim(2));
 end
 disp('Total time:')
 minutes((now - st_0)*24*60)
     
 % Save variables to .nc files 
 disp('Saving')     
 nccreate([fpath '/climatology.nc'],'af_mean','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
 ncwrite( [fpath '/climatology.nc'],'af_mean',af)
 ncwriteatt([fpath '/climatology.nc'],'af_mean','units','none')
 ncwriteatt([fpath '/climatology.nc'],'af_mean','long_name','alkalinity factor monthly clmatology for control fun')
     
 nccreate([fpath '/climatology.nc'],'sf_mean','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
 ncwrite( [fpath '/climatology.nc'],'sf_mean',sf)
 ncwriteatt([fpath '/climatology.nc'],'sf_mean','units','none')
 ncwriteatt([fpath '/climatology.nc'],'sf_mean','long_name','salinity factor monthly clmatology for control fun')


%% Transient Run Alkalinity & Salinity Factors -----------------------------------
 
clearvars -except data_path 
clc

% Find unique years
time  = ncread([data_path '/Ensemble_Members/time.nc'],'time');
yr    = str2num(datestr(time,'yyyy'));
uyr   = unique(yr);
 
% Perform calculations globally in 5 year chunks
yint     = 5;
interval = 12*yint;
yind     = [1:interval:length(time)];  
uyind    = 1:yint:length(uyr);

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim = vinfo.Size;  
 
st_0 = now;
for ii = 1:30
    clearvars -except data_path st_0 ii yr uyr yind uyind interval yint dim time
    
    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now; 

    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    
    % Pre-allocate
    af = NaN(dim);
    sf = NaN(dim);

    % Perform calculations in chunks
    for i = 1:length(yind)
        clear  DIC SSS SST SIO4 PO4 CO2 TA xx A 

        disp(['Start chunk ' num2str(i) ' of ' num2str(length(yind))])   
        st_2 = now;        
        if i < length(yind)
            yint 	 = 5;
            interval = 12*yint;
    	elseif i == length(yind)
            % Change interval for last chunk
            yint = 1;
            interval = length(time) - yind(i) + 1;
    	end

        % Load chunk of data and convert to vector
        SST  = reshape(ncread([fpath  '/sst.nc'],'sst', [1 1 yind(i)], [Inf Inf interval]),[],1);
        SSS  = reshape(ncread([fpath  '/sss.nc'],'sss', [1 1 yind(i)], [Inf Inf interval]),[],1);
        SIO4 = reshape(ncread([fpath '/sio4.nc'],'sio4',[1 1 yind(i)], [Inf Inf interval]),[],1);
        PO4  = reshape(ncread([fpath  '/po4.nc'],'po4', [1 1 yind(i)], [Inf Inf interval]),[],1);
        TA   = reshape(ncread([fpath   '/ta.nc'],'ta',  [1 1 yind(i)], [Inf Inf interval]),[],1);
        CO2  = reshape(ncread([fpath  '/co2.nc'],'co2', [1 1 yind(i)], [Inf Inf interval]),[],1);
        DIC  = reshape(ncread([fpath  '/dic_calc.nc'],'dic_calc', [1 1 yind(i)], [Inf Inf interval]),[],1);
     
        % Define purturbation
        ta_add  = [-0.01 0.01];
        sal_add = [-0.0001 0.0001];
                  
        % Find non-NaN values
        xx = find(~isnan(TA) & ~isnan(DIC) & ~isnan(SIO4));
        % Perform CO2SYS calculations for modern/chemical leverage conditions

        disp('Start CO2SYS calculations')
        % Pre-allocate
        build_af = NaN(length(TA), length(ta_add));
        build_sf = NaN(length(TA), length(ta_add));

        % Calculate sensitivity terms
        for j = 1:length(ta_add)
            % Alkalinity Factor
            A    = CO2SYS(TA(xx) + ta_add(j),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
            dco2 = A(:,19) - CO2(xx);
            build_af(xx,j) = (dco2 ./ CO2(xx)) .* (ta_add(j) ./ TA(xx)).^-1;   
            clear A dco2
         
            % Salinity Factor
            A    = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx) + sal_add(j) ,SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
            dco2 = A(:,19) - CO2(xx);
            build_sf(xx,j) = (dco2 ./ CO2(xx)) .* (sal_add(j) ./ SSS(xx)).^-1; 
            clear A dco2
        end
     
        af(:,:,yind(i):yind(i)+interval-1)  = reshape(nanmean(build_af,2),dim(1),dim(2),[]);
        sf(:,:,yind(i):yind(i)+interval-1)  = reshape(nanmean(build_sf,2),dim(1),dim(2),[]);
     
        disp(['Complete chunk ' num2str(i)])
        minutes((now - st_2)*24*60)
    end  
    disp('Finish Calculations')
    minutes((now - st_1)*24*60)
 
    % Save variables to .nc files 
    disp('Saving')     
 
    nccreate([fpath   '/af.nc'],'af','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/af.nc'],'af',af)
    ncwriteatt([fpath '/af.nc'],'af','units','none')
    ncwriteatt([fpath '/af.nc'],'af','long_name','alkalinity factor')
    
    nccreate([fpath   '/sf.nc'],'sf','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
    ncwrite( [fpath   '/sf.nc'],'sf',sf)
    ncwriteatt([fpath '/sf.nc'],'sf','units','none')
    ncwriteatt([fpath '/sf.nc'],'sf','long_name','salinity factor')
     
    disp(['Complete ensemble member ' num2str(ii)])
    minutes((now - st_1)*24*60) 
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~45 minutes per ensemble member


%% Calculate Biophysical Component w/ and w/o Leverage -----------------------------------
 
clearvars -except data_path 
clc
 
fpath = [data_path '/Ensemble_Members/m'];
% Get final file size from existing variable (ta)
vinfo = ncinfo([fpath '1/ta.nc'],'ta');
dim = vinfo.Size;  

rf_pi    = reshape(repmat(ncread([data_path '/Control_Run/climatology.nc'],'rf_mean'),1,1,1,dim(3)/12),[],1);
af_pi    = reshape(repmat(ncread([data_path '/Control_Run/climatology.nc'],'af_mean'),1,1,1,dim(3)/12),[],1);
clear v M
 
% Cycle through 30 ensemble members
st_0 = now;
for ii = 1%:30   
    st_1 = now;
    disp(['Start ensemble member ' num2str(ii)])
    clearvars -except data_path st_0 st_1 rf_pi af_pi ii dim fpath
    
    % DIC
    dic_calc = ncread([fpath num2str(ii) '/dic_calc.nc'],'dic_calc');
    DIC      = reshape(dic_calc,[],1);
    % Calculate annual mean, repeat annual value over 12 months of each year (yint), convert to vector
    dic      = reshape(repmat(nanmean(reshape(dic_calc,dim(1),dim(2),12,[]),3),1,1,12,1),[],1);
    % Monthly anomalies
    d_dic     = DIC - dic; 
    clear DIC dic_calc
    
    % TA
    ta_in  = ncread([fpath num2str(ii) '/ta.nc'],'ta');
    TA     = reshape(ta_in,[],1);
    % Calculate annual mean, repeat annual value over 12 months of each year (yint), convert to vector
    ta     = reshape(repmat(nanmean(reshape(ta_in,dim(1),dim(2),12,[]),3),1,1,12,1),[],1);
    % Monthly anomalies
    d_ta    = TA - ta; 
    clear TA ta_in
    
    % SSS
    sss_in  = ncread([fpath num2str(ii) '/sss.nc'],'sss');
    SSS     = reshape(sss_in,[],1);
    % Calculate annual mean, repeat annual value over 12 months of each year (yint), convert to vector
    sss     = reshape(repmat(nanmean(reshape(sss_in,dim(1),dim(2),12,[]),3),1,1,12,1),[],1);
    % Monthly anomalies
    d_sal    = SSS - sss; 
    clear SSS sss_in
    
    % CO2 AM values
    CO2_AM = ncread([fpath num2str(ii) '/co2_am.nc'],'co2_am');
    % Repeat annual value over 12 months of each year (yint), convert to vector
    co2_am = reshape(repmat(reshape(CO2_AM,dim(1),dim(2),1,dim(3)/12),1,1,12,1),[],1); 
    
    % Load sensitivity term and convert to vector
    rf_t   = reshape(ncread([fpath num2str(ii) '/rf.nc'],'rf'),[],1);
    af_t   = reshape(ncread([fpath num2str(ii) '/af.nc'],'af'),[],1);  
    sf_t   = reshape(ncread([fpath num2str(ii) '/sf.nc'],'sf'),[],1);    
    clear CO2_AM  

    disp('Start calculations')
    xx  = find(dic > 0 & ta > 0); 
    
    % No-leverage calculations   
    d_co2_nl = NaN(size(d_dic));
    d_co2_nl(xx) = (rf_pi(xx) .* (d_dic(xx) ./ dic(xx)) .* co2_am(xx)) + ...
                   (af_pi(xx) .* (d_ta(xx)  ./  ta(xx)) .* co2_am(xx)) + ...
                   (sf_t(xx)  .* (d_sal(xx) ./ sss(xx)) .* co2_am(xx));
    co2_bp_nl = reshape((d_co2_nl + co2_am),dim(1),dim(2),[]);
      
    % Chemical-leverage calculations   
    d_co2 = NaN(size(d_dic));
    d_co2(xx) = (rf_t(xx) .* (d_dic(xx) ./ dic(xx)) .* co2_am(xx)) + ...
                (af_t(xx) .* (d_ta(xx)  ./  ta(xx)) .* co2_am(xx)) + ...
                (sf_t(xx) .* (d_sal(xx) ./ sss(xx)) .* co2_am(xx));
    co2_bp = reshape((d_co2 + co2_am),dim(1),dim(2),[]);
    
    % Save variables to .nc files 
    disp('Saving')        
    nccreate([fpath num2str(ii) '/co2_bp_nl.nc'],'co2_bp_nl','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5);
    ncwrite([fpath  num2str(ii) '/co2_bp_nl.nc'],'co2_bp_nl',co2_bp_nl);
    ncwriteatt([fpath  num2str(ii) '/co2_bp_nl.nc'],'co2_bp_nl','units','uatm');
    ncwriteatt([fpath  num2str(ii) '/co2_bp_nl.nc'],'co2_bp_nl','long_name','biophysical pco2 component w/o chemical leverage');
    
    nccreate([fpath num2str(ii) '/co2_bp.nc'],'co2_bp','Dimensions',{'lon',dim(1),'lat',dim(2),'time' dim(3)},'DeflateLevel',5);
    ncwrite([fpath  num2str(ii) '/co2_bp.nc'],'co2_bp',co2_bp);
    ncwriteatt([fpath  num2str(ii) '/co2_bp.nc'],'co2_bp','units','uatm');
    ncwriteatt([fpath  num2str(ii) '/co2_bp.nc'],'co2_bp','long_name','biophysical pco2 component w/ chemical leverage');
    clear d_co2_nl d_co2 co2_am xx dic sss ta rf_t af_t sf_f  
    
    disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)
~5 min per ensemble member
 
 
%% Clean Up -----------------------------------
 
clear;close all;clc
