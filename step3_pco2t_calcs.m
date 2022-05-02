%% step3_pco2t_calcs.m
%
% DESCRIPTION:
%	Calculate thermal pCO2 component (pCO2 T) w/ and w/o chemical leverage.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT: 
%	co2_t_nl (microatmosphere) - ensemble member thermal pCO2 component w/o chemical leverage (nl)
%	co2_t    (microatmosphere) - ensemble member thermal pCO2 component w/ chemical leverage
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


%% Calculate Thermal pCO2 Component w/ and w/o Leverage ----------------------------------- 

% Find unique years
time  = ncread([data_path '/Ensemble_Members/time.nc'],'time');
yr    = str2num(datestr(time,'yyyy'));
uyr   = unique(yr); 

% Perform calculations globally in chunks
yint 	 = 10;
interval = 12*yint;
yind     = [1:interval:length(time)]; 
uyind    = 1:yint:length(uyr); 

% Load control run co2_am value, repeat annual value over 12 months of each year (yint), then convert to vector
co2_am_pi_full = reshape(repmat(ncread([data_path '/Control_Run/co2_am.nc'],'co2_am'),1,1,12*yint),[],1);
% Load control run co2_am value, repeat annual value over 12 months for one year, then convert to vector
co2_am_pi_1yr  = reshape(repmat(ncread([data_path '/Control_Run/co2_am.nc'],'co2_am'),1,1,12),[],1);

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim = vinfo.Size;  

st_0 = now;
for ii = 1:30
	clearvars -except data_path st_0 ii time uyr co2_am_pi_1yr co2_am_pi_full yint dim yind uyind

    disp(['Start ensemble member ' num2str(ii)])   
	st_1 = now; 

	fpath = [data_path '/Ensemble_Members/m' num2str(ii)];  

    % Pre-allocate
    co2_t_nl = NaN(dim);
    co2_t    = NaN(dim);

    % Perform calculations in chunks
    for i = 1:length(yind)
    	clear  sst_av sss_av si_av po4_av ta_av co2_am sst co2 A xx var co2_am_pi

    	disp(['Start chunk ' num2str(i) ' of ' num2str(length(yind))])   
    	st_2 = now;

        if i < length(yind)
			yint 	 = 10;
			interval = 12*yint;
			co2_am_pi = co2_am_pi_full;
		elseif i == length(yind)
        	% Change interval for last chunk
        	yint = 1;
        	interval = length(time) - yind(i) + 1;
			co2_am_pi = co2_am_pi_1yr;
        end

	    % Load chunk of data, average data by year, repeat annual mean value over 12 months of each year, then convert to vector
	    sst_av = reshape(repmat(nanmean(reshape(ncread([fpath  '/sst.nc'],'sst', [1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	    sss_av = reshape(repmat(nanmean(reshape(ncread([fpath  '/sss.nc'],'sss', [1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	    si_av  = reshape(repmat(nanmean(reshape(ncread([fpath '/sio4.nc'],'sio4',[1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	    po4_av = reshape(repmat(nanmean(reshape(ncread([fpath  '/po4.nc'],'po4', [1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	    ta_av  = reshape(repmat(nanmean(reshape(ncread([fpath   '/ta.nc'],'ta',  [1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	    co2_av = reshape(repmat(nanmean(reshape(ncread([fpath  '/co2.nc'],'co2', [1 1 yind(i)], [Inf Inf interval]),dim(1),dim(2),12,yint),3),1,1,12,1),[],1);
	   
	    % Load chunk of data, repeat annual mean value over 12 months of each year, then convert to vector
	    co2_am = reshape(repmat(reshape(ncread([fpath '/co2_am.nc'],'co2_am',[1 1 uyind(i)],[Inf Inf yint]), dim(1),dim(2),1,yint),1,1,12,1),[],1);

	    % Load chunk of data, then convert to vector
	    sst    = reshape(ncread([fpath  '/sst.nc'],'sst', [1 1 yind(i)], [Inf Inf interval]),[],1);
	    co2    = reshape(ncread([fpath  '/co2.nc'],'co2', [1 1 yind(i)], [Inf Inf interval]),[],1);
        
	  	% Find non-NaN values
	  	xx = find(ta_av>0 & sss_av>0 & co2_am_pi>0);
	      
	  	% Perform CO2SYS calculations for no-leverage conditions
	  	disp('Start no-leverage CO2SYS calculations')
	    A = CO2SYS(ta_av(xx),co2_am_pi(xx),1,4,sss_av(xx),sst_av(xx),sst(xx),0,0,si_av(xx),po4_av(xx),1,4,1);        
	    var = NaN(size(ta_av));
	    var(xx) = A(:,19) - co2_am_pi(xx) + co2_am(xx);
	    co2_t_nl(:,:,yind(i):yind(i)+interval-1) = reshape(var,dim(1),dim(2),[]);
	
	    clear xx var A 
	
	  	% Find non-NaN values
	  	xx = find(ta_av>0 & sss_av>0 & co2_am>0);
	      
	  	% Perform CO2SYS calculations for chemical-leverage conditions
	  	disp('Start chem-leverage CO2SYS calculations')
	  	A = CO2SYS(ta_av(xx),co2_am(xx),1,4,sss_av(xx),sst_av(xx),sst(xx),0,0,si_av(xx),po4_av(xx),1,4,1);
	    var = NaN(size(ta_av));
	    var(xx) = A(:,19);
	    co2_t(:,:,yind(i):yind(i)+interval-1) = reshape(var,dim(1),dim(2),[]);

		disp(['Complete chunk ' num2str(i)])
		minutes((now - st_2)*24*60)
	end
    disp('Finish Calculations')
    minutes((now - st_1)*24*60)

	% Save variables to .nc files 
	disp('Saving') 

	nccreate([fpath   '/co2_t_nl.nc'],'co2_t_nl','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
	ncwrite([fpath    '/co2_t_nl.nc'],'co2_t_nl',co2_t_nl)
	ncwriteatt([fpath '/co2_t_nl.nc'],'co2_t_nl','units','uatm')
	ncwriteatt([fpath '/co2_t_nl.nc'],'co2_t_nl','long_name','thermal pco2 component w/o chemical leverage')

	nccreate([fpath   '/co2_t.nc'],'co2_t','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)},'DeflateLevel',5)
	ncwrite([fpath    '/co2_t.nc'],'co2_t',co2_t)
	ncwriteatt([fpath '/co2_t.nc'],'co2_t','units','uatm')
	ncwriteatt([fpath '/co2_t.nc'],'co2_t','long_name','thermal pco2 component w/ chemical leverage') 

	disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
	minutes((now - st_1)*24*60)	
end
disp('Total time:')
minutes((now - st_0)*24*60)
% ~13 hours


%% Clean Up -----------------------------------

clear;close all;clc
