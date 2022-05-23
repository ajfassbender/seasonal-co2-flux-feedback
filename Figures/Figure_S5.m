%% Figure S5
%
% DESCRIPTION:
%   Model output vs. reconstructed pCO2 (prior to bias correction) x-y plot 
%   Ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   cmocean: https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
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

time = ncread([data_path '/Ensemble_Members/time.nc'],'time'); 
lat  = ncread([data_path '/Ensemble_Members/lat.nc'],'lat'); 
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 

fpath    = [data_path '/Ensemble_Members/m1/'];
co2      = ncread([fpath 'co2.nc'], 'co2');
co2_totl = ncread([fpath 'co2_totl.nc'], 'co2_totl');


%% Make Figure -----------------------------------

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 10 6.5 2])
set(gcf,'Resize','off')

subplot(121)
x = reshape(co2,[],1);
y = reshape(co2_totl,[],1);
z = find(isfinite(x)==1 & isfinite(y)==1);

% [r,p] = corrcoef(x(z),y(z))
% r = 0.9997
% p < 0.01

h = histogram2(x(z),y(z),'NumBins',100,'DisplayStyle','tile','ShowEmptyBins','off');
hold on
plot([-200 1500],[-200 1500],'r'); hold on
grid on;h = colorbar;
col = cmocean('thermal',100);
set(gca,'colormap',col)
set(gca,'ColorScale','log')
xlim([-200 1500]);ylim([-200 1500])
xlabel('Model Output');
ylabel('Reconstructed')
title('\itp\rmCO_2 (\muatm)')
ylabel(h,'histogram frequency')
caxis([1 10^6])
set(h,'XTick',[10^0 10^2 10^4 10^6],'xticklabel',{'10^0','10^2','10^4','10^6'})

clear x y z col h n 

subplot(122) 
lat_ind = fliplr(0.5:10:90.5);
col     = cmocean('thermal',9);
m       = 1;
for i = 1:length(lat_ind)-1;
	clear x y z ii
	z = find(lat_ind(i)>=abs(lat) & abs(lat)>lat_ind(i+1));
	x = reshape(co2(:,z,:),[],1);
	if nansum(x(:))~=0
		y  = reshape(co2_totl(:,z,:),[],1);
		ii = find(isfinite(x)==1 & isfinite(y)==1);
		plot(x(ii),y(ii),'.','color',col(m,:),'markersize',.1);hold on
		m = m+1;
		clear x y ii z
	end
	disp(['Done with latitude band ' num2str(i) ' of ' num2str(length(lat_ind)-1)])
end
clear x y ii co2 co2_totl
plot([-200 1500],[-200 1500],'r'); 
xlim([-200 1500]);ylim([-200 1500])
xlabel('Model Output');
ylabel('Reconstructed')
title('\itp\rmCO_2 (\muatm)')

cind = fliplr({'  0.5-10.5','10.5-20.5','20.5-30.5','30.5-40.5','40.5-50.5','50.5-60.5','60.5-70.5','70.5-80.5','80.5-90.5'});
inc  = 0+(1/18):1/9:1;
h    = colorbar;
set(h,'XTick',inc,'xticklabel',cind)
ylabel(h,'|Latitude|')
grid on;set(gca,'colormap',col)


%% Subplot Letters -----------------------------------

let = {'(a)','(b)'};
for n = 1:length(let)
    subplot(1,2,n) ;hold on    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.05*xr,yylim(2)-0.15.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S5.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S5'])


%% Clean Up -----------------------------------

clear;close all;clc
