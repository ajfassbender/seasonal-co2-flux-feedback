%% Figure 9
%
% DESCRIPTION:
%   Flux bias from using annual mean CO2
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   subtightplot: https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
%   rgb: https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2/
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


%% Set Regional Bounaries and Compute Area Weighted Averages -----------------------------------

% REGION = [Lon_min Lon_max Lat_min Lat_max]
 v_ind = {'na','np','st','eq','so'};
 r_ind = [ 50 60 315 350 ;
           42 55 150 185 ;
           30 40 145 225 ;
          -10 10 210 275 ;
          -60 -40 1 360 ];

fpath = [data_path '/Ensemble_Members'];
lat   = ncread([fpath '/lat.nc'],  'lat');   
lon   = ncread([fpath '/lon.nc'],  'lon');  
time  = ncread([fpath '/time.nc'], 'time'); 
area  = transpose(getArea(-89.5:1:89.5,.5:359.5));

flux_co2_av = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_annual_mean_co2');
flux_co2    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_l');

for j = 1:length(v_ind)
    x1 = find(r_ind(j,1) <= lat & lat <= r_ind(j,2));
    x2 = find(r_ind(j,3) <= lon & lon <= r_ind(j,4));
    flux_avg(:,j+1) = nansum(nansum(flux_co2_av(x2,x1) .*  area(x2,x1))) .* 12 ./ (10^15);% Pg C
    flux(:,j+1)     = nansum(nansum(flux_co2(x2,x1)    .*  area(x2,x1))) .* 12 ./ (10^15);% Pg C
end

% Add Total Uptake
flux_avg(:,1)    = nansum(nansum(flux_co2_av .*  area)) .* 12 ./ (10^15);% Pg C
flux(:,1)        = nansum(nansum(flux_co2    .*  area)) .* 12 ./ (10^15);% Pg C


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

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];
x = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2 = [transpose(llon) .5];

clear N XX
N = transpose(flux_co2 - flux_co2_av);
XX = [N N(:,1)];

figure(2);clf
set(gcf,'Units','Inches');
set(gcf,'position',[0 12 6.5 2])
set(gcf,'Resize','off')

subtightplot(1,2,1,0,0,[0 0])

for l=1:6
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
        'ytick',[],'yticklabels',[],'linest','-','color',[.2 .2 .2]);
    if l<3
        h=m_pcolor(LON2,lat,XX);hold on;%shading interp
        set(h, 'EdgeColor','none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.75 .75 .75],'xticklabels',[],'yticklabels',[]);
    else
        h=m_pcolor(lon,lat,N);hold on;%shading interp
        set(h, 'EdgeColor', 'none');
        m_grid('tickdir','in','linestyle','none','backcolor',[.75 .75 .75],'xticklabels',[],'yticklabels',[]);
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
T= title({'Impact of Using Monthly Data','(mol C m^{-2})'});set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
h = colorbar;grid on
cc = cmocean('balance',9);
colormap(gca,cc);caxis([-150 150])
set(h,'position',[0.0743 0.2499 0.0176 0.2333])
fig=gcf;
set(findall(fig,'-property','Fontname'),'Fontname','Times','fontsize',9,'fontweight','normal')%% Plot


%% Cumulative Flux Differences -----------------------------------

% Area Weighted Percentage Differences
tot = 100.*(flux-flux_avg)./abs(flux);
round(tot) % -14    17    28   -36    -2     1

% Error propagagion 
flux_avf_sd = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_annual_mean_co2_stdev');
flux_sd = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_l_sd');
v = (flux_avf_sd.^2 + flux_sd.^2).^.5 % Pg C
diffy = flux(1)-flux_avg(1);
tot(1) .* ((v./(diffy)).^2 + (flux_sd./flux(1)).^2 ) .^.5 % 0.37 percent
Pg_vals = num2str(transpose(round(flux - flux_avg)));

figure(2);hold on
subtightplot(1,2,2,.08,[.2 0.225],[.12 .1])
col = [rgb('black');rgb('blueviolet');rgb('dodgerblue'); rgb('orange');rgb('deeppink');rgb('limegreen')];
H = bar(1:6,tot,'FaceColor','flat','BarWidth',.8,'EdgeColor','none');hold on
H.CData = rgb('black');
ylim([-50 50]);ylabel('%')
title({'Cumulative Uptake Difference','\itp\rmCO_{2 Month} - \itp\rmCO_{2 Avg.}'})
set(gca,'xticklabel',{'Global','SP NA','SP NP','ST NP','EqPac','SO'},'Xticklabelrotation',40)

r  = [-1 1 1 -1 -1 -1];
text(.6:6,round(tot+7.*r),[char(Pg_vals)],'color','k');
set(gca,'ytick',-50:25:50)


%% Subplot Letters -----------------------------------

subtightplot(1,2,1,0,0,[0 0])
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)+0.075.*xr,yylim(2)+0.1.*yr,'(a)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',9,'fontweight','normal','fontname','Arial')

subtightplot(1,2,2,.08,[.2 0.225],[.12 .1])
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)-0.22.*xr,yylim(2)+0.1.*yr,'(b)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',9,'fontweight','normal','fontname','Arial')


%% Save -----------------------------------

fig=gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_9.pdf'],'Resolution',300)
savefig([fig_path '/Figure_9'])


%% Clean up -----------------------------------

clear;close all;clc
