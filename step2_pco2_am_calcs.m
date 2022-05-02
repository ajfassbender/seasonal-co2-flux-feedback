%% step2_pco2_am_calcs.m
% 
% DESCRIPTION:
%    Calculate pCO2_AM for the control run and ensemble members.
%
% INPUT:
%   data_path: directory containing reprocessed .nc files
%
% OUTPUT: 
%    co2_am (microatmosphere) - for control run and ensemble members
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

%% Control Run pCO2_AM Calculations -----------------------------------

fpath = [data_path '/Control_Run/'];

vinfo = ncinfo([fpath 'ta.nc'],'ta');
dim = vinfo.Size;  

% Calculate annual mean values from climatologies and reshape to vectors
DIC  = reshape(nanmean(ncread([fpath 'climatology.nc'],'dic_calc_mean'),3),[],1);
TA   = reshape(nanmean(ncread([fpath 'climatology.nc'],'ta_mean'),3),[],1);
SSS  = reshape(nanmean(ncread([fpath 'climatology.nc'],'sss_mean'),3),[],1);
SST  = reshape(nanmean(ncread([fpath 'climatology.nc'],'sst_mean'),3),[],1);
SIO4 = reshape(nanmean(ncread([fpath 'climatology.nc'],'sio4_mean'),3),[],1);
PO4  = reshape(nanmean(ncread([fpath 'climatology.nc'],'po4_mean'),3),[],1);

% Find non-NaN values
xx = find(~isnan(TA) & ~isnan(SSS) & ~isnan(SIO4));

% Start CO2SYS calculations
A    = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
var  = NaN(size(DIC));
Var(xx,1) = A(:,19);
co2_am    = reshape(Var,dim(1),dim(2));

% Save variable to .nc file 
disp('Saving')     
nccreate([fpath 'co2_am.nc'],'co2_am','Dimensions',{'lon',dim(1),'lat',dim(2)},'DeflateLevel',5)
ncwrite( [fpath 'co2_am.nc'],'co2_am',co2_am)
ncwriteatt([fpath 'co2_am.nc'],'co2_am','units','uatm')
ncwriteatt([fpath 'co2_am.nc'],'co2_am','long_name','annual mean pco2 computed from annual mean values of dic, ta, sss, sst, po4, and sio4')


%% Transient Run pCO2_AM Calculations -----------------------------------

clearvars -except data_path

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim = vinfo.Size; clear vinfo

st_0 = now;
for ii = 1:30
	clearvars -except data_path st_0 ii dim

	st_1 = now;

	fpath = [data_path '/Ensemble_Members/m' num2str(ii)];
    disp(['Start ensemble member ' num2str(ii)])

	% Calculate annual means
	dic  = nanmean(reshape(ncread([fpath '/dic_calc.nc'],'dic_calc'),dim(1),dim(2),12,dim(3)/12),3);
	ta   = nanmean(reshape(ncread([fpath '/ta.nc'],'ta'),dim(1),dim(2),12,dim(3)/12),3);
	sss  = nanmean(reshape(ncread([fpath '/sss.nc'],'sss'),dim(1),dim(2),12,dim(3)/12),3);
	sst  = nanmean(reshape(ncread([fpath '/sst.nc'],'sst'),dim(1),dim(2),12,dim(3)/12),3);
	sio4 = nanmean(reshape(ncread([fpath '/sio4.nc'],'sio4'),dim(1),dim(2),12,dim(3)/12),3);
	po4  = nanmean(reshape(ncread([fpath '/po4.nc'],'po4'),dim(1),dim(2),12,dim(3)/12),3);

    % Reshape to vectors
    DIC  = reshape(dic, [],1);  clear dic
    TA   = reshape(ta,  [],1);  clear ta
    SSS  = reshape(sss, [],1);  clear sss
    SST  = reshape(sst, [],1);  clear sst
    SIO4 = reshape(sio4,[],1);  clear sio4
    PO4  = reshape(po4, [],1);  clear po4
    
	% Find non-NaN values
	xx   = find(~isnan(TA) & ~isnan(SSS) & ~isnan(SIO4));

    disp('Start CO2SYS calculations')
	A    = CO2SYS(TA(xx),DIC(xx),1,2,SSS(xx),SST(xx),SST(xx),0,0,SIO4(xx),PO4(xx),1,4,1);
	var  = NaN(size(DIC));
	Var(xx,1) = A(:,19); clear A xx
	co2_am    = reshape(Var,dim(1),dim(2),dim(3)/12);
        
    % Save variable to .nc file 
    disp('Saving') 
	nccreate([fpath   '/co2_am.nc'],'co2_am','Dimensions',{'lon',dim(1),'lat',dim(2),'time',dim(3)/12},'DeflateLevel',5)
	ncwrite( [fpath   '/co2_am.nc'],'co2_am',co2_am)
	ncwriteatt([fpath '/co2_am.nc'],'co2_am','units','uatm')
	ncwriteatt([fpath '/co2_am.nc'],'co2_am','long_name','annual mean pco2 computed from annual mean values of dic, ta, sss, sst, po4, and sio4')
        
    disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
    minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Clean Up -----------------------------------

clear;close all;clc

