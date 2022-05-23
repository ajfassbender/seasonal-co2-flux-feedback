%% Figure 6
%
% DESCRIPTION:
%   Impact of chemical leverage on seasonal pCO2 values and fluxes &
%   local seasonal wind speed vs. dpCO2 correlations 
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
%   tight_subplot: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
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

%% Load Variables -----------------------------------

fpath = [data_path '/Ensemble_Members/'];
vars  = {'diff_co2_summer','diff_co2_winter','diff_flux_summer','diff_flux_winter','diff_aflux'};
for i = 1:length(vars)
    uvar = char(vars(i));
    x = ncread([fpath 'Figure6_impact_of_chem_lev.nc'],uvar);
    eval([uvar '= x;']);
end

vars = {'R','Pval','R_nl','Pval_nl'};
for i = 1:length(vars)
    uvar = char(vars(i));
    x = ncread([fpath 'Figure6_wind_dco2_correlations.nc'],uvar);
    eval([uvar '= x;']);
end


%% Define Regional Boxes (to align w/ mmap projection) -----------------------------------

% REGION = [Lon_min Lon_max Lat_min Lat_max];
% NASP   = [315 350 140 155];
% NPSP   = [150 185 132 145];
% NPST   = [145 225 120 130];
% Eq     = [210 275  80 100];
% SO     = [  1 360  40  60];

X = NaN(5,5,6);
Y = NaN(5,5,6);

% North Atlantic SubPolar Box
Y(1,1:5,1) = [50 50 65 65 50];
X(1,1:5,1) = [315-360 350-360 350-360 315-360 315-360];
   
% Southern Ocean Lines
Y(1,1:2,2) = [-40 -40];  
X(1,1:2,2) = [-65 25.5];
Y(2,1:2,2) = [-60 -60];  
X(2,1:2,2) = [-75 25.5];
 
% Southern Ocean Lines    
Y(1,1:2,3) = [-40 -40];
X(1,1:2,3) = [ 24 145];
Y(2,1:2,3) = [-60 -60];
X(2,1:2,3) = [ 24 145];

% Southern Ocean Lines
Y(1,1:2,5) = [-40 -40];
X(1,1:2,5) = [142 287];
Y(2,1:2,5) = [-60 -60];  
X(2,1:2,5) = [142 295];

% Equatorial Pacific Box (do not draw line at equator)
Y(3,1:2,5) = [-10 -10];  
X(3,1:2,5) = [210 275];
Y(4,1:2,5) = [-10 0];  
X(4,1:2,5) = [210 210];
Y(5,1:2,5) = [-10 0];  
X(5,1:2,5) = [275 275];

% North Pacific Boxes
Y(1,1:5,6) = [42 42 55 55 42];
X(1,1:5,6) = [150 185 185 150 150];
Y(2,1:5,6) = [30 30 40 40 30];
X(2,1:5,6) = [145 225 225 145 145];

% Equatorial Pacific Box (do not draw line at equator)
Y(3,1:2,6) = [10 10];  
X(3,1:2,6) = [210 275];
Y(4,1:2,6) = [10 0];  
X(4,1:2,6) = [210 210];
Y(5,1:2,6) = [10 0];  
X(5,1:2,6) = [275 275];


%% Make Figure -----------------------------------

lat   = ncread([fpath '/lat.nc'],  'lat');   
lon   = ncread([fpath '/lon.nc'],  'lon'); 

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[1 10 6.5 5])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(3,2,[0.06 0],[0 .08],[0 0]);

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

x = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2 = [transpose(llon) .5];


%% Seasonal Max & Min Delta pCO2 -----------------------------------

s_var =  transpose(diff_co2_summer);
Sum = [s_var s_var(:,1)];

w_var =  transpose(diff_co2_winter);
Wint = [w_var w_var(:,1)];

axes(ha(1));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,Sum);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,s_var);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end          
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Impact of Chem. Lev. on','2090s Mean Summer \itp\rmCO_2 (\muatm)'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.073 0.6794 0.0175 0.1])
bn = -20:10:80;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn),'pivot',0);
colormap(gca,cc);


axes(ha(2));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,Wint);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,w_var);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end   
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Impact of Chem. Lev. on','2090s Mean Winter \itp\rmCO_2 (\muatm)'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.572 0.6794 0.0175 0.1])
bn = -80:10:20;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn),'pivot',0);
colormap(gca,cc);


%% Max & Min Flux Anomalies -----------------------------------

mx_f = transpose(diff_flux_summer); %mol/m2/sesaon
mx_F = [mx_f mx_f(:,1)];

mn_f = transpose(diff_flux_winter);
mn_F = [mn_f mn_f(:,1)];

axes(ha(3));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,mx_F);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,mx_f);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end   
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Impact of Chem. Lev. on','2090s \SigmaFlux_{Summer} (mol m^{-2})'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.073 0.353 0.0175 0.1])
bn = -2:0.4:2;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn));
colormap(gca,cc);


axes(ha(4));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,mn_F);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,mn_f);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end   
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Impact of Chem. Lev. on','2090s Winter \SigmaFlux_{Winter} (mol m^{-2})'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.572 0.353 0.0175 0.1])
bn = -2:0.4:2;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn));
colormap(gca,cc);


%% Annual Cumulative Flux Difference -----------------------------------

daflux = transpose(diff_aflux);
Daflux = [daflux daflux(:,1)];

axes(ha(5));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,Daflux);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,daflux);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end   
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Impact of Chem. Lev. on','2090s \SigmaFlux_{Annual} (mol m^{-2} yr^{-1})'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.073 0.0268 0.0175 0.1])
bn = -1.5:.3:1.5;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn),'pivot',0);
colormap(gca,cc);


%% Wind Speed & Correlation -----------------------------------

puse = transpose(nanmean(Pval,3));
mask = puse>0.05;

ruse = transpose(nanmean(R,3));
Ruse = [ruse ruse(:,1)];

axes(ha(6));hold on
for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.4 .4 .4]);
    if l<3
        h=m_pcolor(LON2,lat,Ruse);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        clear xx yy latuse lonuse Xd Yd
        [xx,yy] = meshgrid(llon,lat);
        latuse = reshape(yy,[],1);
        lonuse = reshape(xx,[],1);
        [Xd,Yd] = m_ll2xy(lonuse,latuse);
        stipple(reshape(Xd,180,360), reshape(Yd,180,360),mask,'color','k','density',90,'markersize',2)
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,ruse);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        clear xx yy latuse lonuse Xd Yd
        [xx,yy] = meshgrid(lon,lat);
        latuse = reshape(yy,[],1);
        lonuse = reshape(xx,[],1);
        [Xd,Yd] = m_ll2xy(lonuse,latuse);
        stipple(reshape(Xd,180,360), reshape(Yd,180,360),mask,'color','k','density',90,'markersize',2)
        m_grid('tickdir','in','linestyle','none','backcolor',[.5 .5 .5],'xticklabels',[],'yticklabels',[]);
    end

    for j = 1:5
        y = Y(j,:,l);
        x = X(j,:,l);
        z = isfinite(x)==1;
        if isempty(z) == 1
        else
            [xx,yy] = m_ll2xy(x(z),y(z));
            plot(xx, yy, 'k', 'linewidth',1.5)
        end     
    end   
    m_coast('patch',[.5 .5 .5]);
end
ax = gca; 
ax.SortMethod = 'childorder'; 

set(gca,'xlimmode','auto','ylimmode','auto');
T= title({'Seasonal Correlation Coefficient','2090s \Delta\itp\rmCO_{2 Sea-Air} vs. Wind Speed'},'fontweight','normal');
set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;
set(h,'position',[0.572 0.0268 0.0175 0.1])
bn = -1:.2:1;
caxis([min(bn) max(bn)]);
cc=cmocean('balance',length(bn),'pivot',0);
colormap(gca,cc);


%% Subplot letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for n = 1:6
    axes(ha(n));hold on
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.1.*xr,yylim(2)-.15.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',11,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig=gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_6.pdf'],'Resolution',300)
savefig([fig_path '/Figure_6'])


%% Clean up -----------------------------------

clear;close all;clc
