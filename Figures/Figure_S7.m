%% Figure S7
%
% DESCRIPTION:
%   Cumulative uptake difference over time between model output and offline flux calculations.
%   Ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   rgb: https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2/
%
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;
fig_path  = cd;


%% Load Data -----------------------------------

area   = transpose(getArea(-89.5:1:89.5,.5:359.5));
time   = ncread([data_path '/Ensemble_Members/time.nc'], 'time');   

fpath    = [data_path '/Ensemble_Members/m1/'];
flux_l   = ncread([fpath 'flux_keff.nc'],'flux_keff'); 
flux_nl  = ncread([fpath 'flux_keff_nl.nc'],'flux_keff_nl'); 
mflux    = ncread([fpath 'mflux.nc'],'mflux'); 
flux_av  = ncread([fpath 'flux_keff_annual_mean_co2.nc'],'flux_keff_annual_mean_co2'); 
flux_offline   = ncread([data_path '/Ensemble_Members/flux_offline_EM1.nc'],'flux_offline'); 

% Ice Mask has not been applied
mflux_mo        = squeeze(nansum(nansum(-mflux  .* area)) .* 365/12 .* 12 ./ 10^15);% Pg C
flux_l_mo       = squeeze(nansum(nansum(flux_l  .* area)) .* 365/12 .* 12 ./ 10^15);% Pg C
flux_nl_mo      = squeeze(nansum(nansum(flux_nl .* area)) .* 365/12 .* 12 ./ 10^15);% Pg C
flux_av_mo      = squeeze(nansum(nansum(flux_av .* area)) .* 365/12 .* 12 ./ 10^15);% Pg C
flux_offline_mo = squeeze(nansum(nansum(flux_offline .* area)) .* 365/12 .* 12 ./ 10^15);% Pg C


%% Make Figure -----------------------------------

figure(2);clf
set(gcf,'Units','Inches');
set(gcf,'position',[1 6 4 2.5])
set(gcf,'Resize','off')

subplot(211)
plot(time,smooth(flux_av_mo,120,'lowess'),'k','color',rgb('black'),'linewidth',1.5);hold on
plot(time,smooth(mflux_mo,120,'lowess'),'k','color',rgb('dodgerblue'),'linewidth',1.5);hold on
plot(time,smooth(flux_offline_mo,120,'lowess'),'k','color',rgb('gray'),'linewidth',1.5);
plot(time,smooth(flux_l_mo,120,'lowess'),'--k','color',rgb('orange'),'linewidth',1.5);
plot(time,smooth(flux_nl_mo,120,'lowess'),'k','color',rgb('blueviolet'),'linewidth',1.5);
grid on;datetick('x','yyyy')
xlim([time(1) time(end)])
ylabel('PgC mo^{-1}');ylim([-.6 0])
set(gca,'ytick',-.6:.2:0)
title('10-Year Smoothed Carbon Uptake')

text(time(830),-0.5,'Flux_{Avg.}','color',rgb('black'))
text(time(30), -0.5,'Flux_{output}','color',rgb('dodgerblue'))
text(time(430),-0.5,'Flux_{offline}','color',rgb('gray'))
text(time(30), -0.3,'Flux_{k_{eff}}','color',rgb('orange'))
text(time(430),-0.3,'Flux_{\0 k_{eff}}','color',rgb('blueviolet'))

subplot(212)
plot(time,smooth(cumsum(flux_av_mo),120,'lowess'),'k','color',rgb('black'),'linewidth',1.5);hold on
plot(time,smooth(cumsum(mflux_mo),120,'lowess'),'k','color',rgb('dodgerblue'),'linewidth',1.5);hold on
plot(time,smooth(cumsum(flux_offline_mo),120,'lowess'),'k','color',rgb('gray'),'linewidth',1.5);
plot(time,smooth(cumsum(flux_l_mo),120,'lowess'),'--k','color',rgb('orange'),'linewidth',1.5);
plot(time,smooth(cumsum(flux_nl_mo),120,'lowess'),'k','color',rgb('blueviolet'),'linewidth',1.5);
grid on;datetick('x','yyyy')
xlim([time(1) time(end)])
ylabel('PgC')
title('Cumulative Carbon Uptake')


%% Subplot Letters -----------------------------------

let = {'(a)','(b)'};
for n = 1:length(let)
    subplot(2,1,n);hold on    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(2)-0.075*xr,yylim(2)-0.225.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S7.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S7'])


%% Clean Up -----------------------------------

clear;close all;clc
