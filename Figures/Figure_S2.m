%% Figure S2 
%
% DESCRIPTION:
%   Difference bewtween annually averaged pCO2 values and pCO2_AM values for the 1950s and 2090s.
%   Use ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
%   tight_subplot: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
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
 
lat  = ncread([data_path '/Ensemble_Members/lat.nc'],  'lat');   
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],  'lon');
time = ncread([data_path '/Ensemble_Members/time.nc'],'time');

yr   = str2num(datestr(time,'yyyy'));
uyr  = unique(yr);
uyr_50   = find(1950 <= uyr & uyr < 1960);
uyr_90   = find(1990 <= uyr & uyr < 2100);

fpath = [data_path '/Ensemble_Members/m1/'];
co2_am_50s = nanmean(ncread([fpath  'co2_am.nc'], 'co2_am', [1 1 uyr_50(1)],[Inf Inf length(uyr_50)]),3);
co2_am_90s = nanmean(ncread([fpath  'co2_am.nc'], 'co2_am', [1 1 uyr_90(1)],[Inf Inf length(uyr_90)]),3);

dim = size(co2_am_90s);
co2_avg    = squeeze(nanmean(reshape(ncread([fpath  'co2.nc'], 'co2'), dim(1),dim(2),12,[]),3));
co2_50s = nanmean(co2_avg(:,:,uyr_50),3);
co2_90s = nanmean(co2_avg(:,:,uyr_90),3);

%%  1950s and 2090s difference in CO2 Av vs. CO2 AM
dco2_50s = co2_50s - co2_am_50s;
dco2_90s = co2_90s - co2_am_90s;


%% Make Figure -----------------------------------

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 12 6.5 1.75])
clear ha;ha = tight_subplot(1,3,[0 0],[0.05 .05],[0 0]);

%map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

var = {'dco2_50s', 'dco2_90s'};
crng1 = [-5 -5];
crng2 = [5 5];
cint  = [1 1];
tvar  = {'1950s: \itp\rmCO_{2 Avg.} - \itp\rmCO_{2 AM}','2090s: \itp\rmCO_{2 Avg.} - \itp\rmCO_{2 AM}'};

x    = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

for n = 1:length(var)
  
    clear xvar
    eval(['xvar = ' char(var(n)) ';']);
    xvar = transpose(xvar);
    clear X;X = [xvar xvar(:,1)];
    
    axes(ha(n));hold on
    
    for l=1:6
        m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
        m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
            'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
        if l<3
            h=m_pcolor(LON2,lat,X);hold on;%shading interp
            set(h, 'EdgeColor', 'none');
            m_coast('patch',[.5 .5 .5]);
            m_grid('tickdir','in','linestyle','none','backcolor',[.75 .75 .75],'xticklabels',[],'yticklabels',[]);
        else
            h=m_pcolor(lon,lat,xvar);hold on;%shading interp
            set(h, 'EdgeColor', 'none');
            m_coast('patch',[.5 .5 .5]);
            m_grid('tickdir','in','linestyle','none','backcolor',[.75 .75 .75],'xticklabels',[],'yticklabels',[]);
        end
    end
    
    set(gca,'xlimmode','auto','ylimmode','auto');
    T= title(tvar(n));set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h = colorbar('southoutside');
    
    bn = crng1(n):cint(n):crng2(n);
    caxis([min(bn) max(bn)]);cc = cmocean('balance',length(bn));
    colormap(gca,cc);
    
    if n == 1
        set(h,'position',[0.0877 0.2014 0.1360 0.0470])
    elseif n == 2
        set(h,'position',[0.43 0.2014 0.1360 0.0470])
    elseif n == 3
        set(h,'position',[0.7649 0.2014 0.1360 0.0470])
    end
end


%% Subplot letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for n = 1:2
    axes(ha(n));hold on
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.05*xr,yylim(2)-0.125.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S2.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S2'])


%% Clean up -----------------------------------

clear;close all;clc
