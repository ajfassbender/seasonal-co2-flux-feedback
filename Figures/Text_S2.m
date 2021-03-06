%% Text S2
%
% DESCRIPTION:
%   Ensemble mean correlation coefficient for pCO2 output vs. recalculated pCO2 (w and w/o bias-correction) 
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
%   area_weighted_mean: computes area weighted mean
%   area_weighted_std: computes area weighted standard deviation
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;

%% pCO2 correlations -----------------------------------

area = transpose(getArea(-89.5:1:89.5,.5:359.5));

st_0 = now;
for ii = 1:30
    clearvars -except data_path fpath ii st_0 area r_coeff av_bias std_bias r_coeff2 av_bias2 std_bias2
    
    fpath = [data_path '/Ensemble_Members/m' num2str(ii)];

    disp(['Start ensemble member ' num2str(ii)])   
    st_1 = now;

    co2      = ncread([fpath  '/co2.nc'],'co2'); 
    co2_corr = ncread([fpath  '/co2_corr.nc'],'co2_corr');

	% Correlations w/ bias correction
	x = reshape(co2,[],1);
	y = reshape(co2_corr,[],1);
	z = find(isfinite(x)==1 & isfinite(y)==1);
	[r,p] = corrcoef(x(z),y(z));
	r_coeff(ii,1) = r(1,2);
	
	diffy = co2 - co2_corr; clear M
	av_bias(ii,1)  = area_weighted_mean(diffy,-89.5:1:89.5,.5:359.5);
	std_bias(ii,1) = area_weighted_std( diffy,-89.5:1:89.5,.5:359.5);
	clear co2_corr x y z r p diffy 

	% Correlations w/o bias correction
	co2_totl = ncread([fpath  '/co2_totl.nc'],'co2_totl');

	x = reshape(co2,[],1);
	y = reshape(co2_totl,[],1);
	z = find(isfinite(x)==1 & isfinite(y)==1);
	[r,p] = corrcoef(x(z),y(z));
	r_coeff2(ii,1) = r(1,2);
	
	diffy = co2 - co2_totl; clear M
	av_bias2(ii,1)  = area_weighted_mean(diffy,-89.5:1:89.5,.5:359.5);
	std_bias2(ii,1) = area_weighted_std( diffy,-89.5:1:89.5,.5:359.5);

	disp(['Complete ensemble member ' num2str(ii)])
	minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)

% Correlations w/ bias correction
nanmean(r_coeff)  %   0.9998
nanmean(av_bias)  %   2.4118e-04
nanmean(std_bias) %   2.4762

% Correlations w/o bias correction
nanmean(r_coeff2)  %    0.9997
nanmean(av_bias2)  %   -1.8232
nanmean(std_bias2) %    3.0349


%% CO2 Flux correlations  -----------------------------------

area = transpose(getArea(-89.5:1:89.5,.5:359.5));

st_0 = now;
for ii = 1:30
	clearvars -except data_path fpath ii st_0 area r_coeff av_bias std_bias r_coeff2 av_bias2 std_bias2

	fpath = [data_path '/Ensemble_Members/m' num2str(ii)];

	disp(['Start ensemble member ' num2str(ii)])   
	st_1 = now;

	flux_keff = ncread([fpath  '/flux_keff.nc'],'flux_keff'); 
	mflux     = -1.*ncread([fpath  '/mflux.nc'],'mflux');

	% Correlations w/ bias correction
	x = reshape(mflux,[],1);
	y = reshape(flux_keff,[],1);
	z = find(isfinite(x)==1 & isfinite(y)==1);
	[r,p] = corrcoef(x(z),y(z));
	r_coeff(ii,1) = r(1,2);

	disp(['Complete ensemble member ' num2str(ii) ' CO2Sys Calcs'])
	minutes((now - st_1)*24*60)
end
disp('Total time:')
minutes((now - st_0)*24*60)

% Correlations
nanmean(r_coeff)  % 0.9872

%% Clean up -----------------------------------

% clear;close all;clc
