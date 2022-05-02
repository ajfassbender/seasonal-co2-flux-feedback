%% Figure S3
%
% DESCRIPTION:
%   pCO2 AM, RF, and AF changes from the pre-industrial to the 2090s.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
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

fpath = [data_path '/Ensemble_Members'];
vars  = {'rf_90s','af_90s','co2_am_90s','rf_90s_sd','af_90s_sd','co2_am_90s_sd'};
for i = 1:length(vars)
    uvar = char(vars(i));
    x = ncread([fpath '/FigureS3_2090s_Ensemble_Mean_RF_AF_CO2_AM.nc'],uvar);
    eval([uvar ' = x;']);
end

% Control Run Data
fpath = [data_path '/Control_Run/'];
rf_pi = nanmean(ncread([fpath 'climatology.nc'], 'rf_mean'),3);
af_pi = nanmean(ncread([fpath 'climatology.nc'], 'af_mean'),3);
co2_am_pi = nanmean(ncread([fpath 'co2_am.nc'], 'co2_am'),3);


%% Make Figure -----------------------------------

lon = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 
lat = ncread([data_path '/Ensemble_Members/lat.nc'],'lat'); 

dco2 = co2_am_90s - co2_am_pi;
drf  = rf_90s - rf_pi;
daf  = af_90s - af_pi;

% How much larger is the longterm change than seasonal variations?
area = transpose(getArea(-89.5:1:89.5,0.5:1:359.5));
marea = nanmean(nanmean(area));
[nanmean(nanmean(dco2.*area))./marea nanmean(nanmean(co2_am_90s_sd.*area))./marea]
[nanmean(nanmean(daf .*area))./marea nanmean(nanmean(af_90s_sd.*area))./marea]
[nanmean(nanmean(drf .*area))./marea nanmean(nanmean(rf_90s_sd.*area))./marea]


var = {'co2_am_90s','rf_90s','af_90s','co2_am_90s_sd','rf_90s_sd','af_90s_sd','dco2','drf','daf'};
crng1 = [800  10  -20   25   0.00    0.00     525     2.0   -8.0 ];
crng2 = [900  20  -10   35   1.20    1.20     625     8.0   -2.0 ];
cint  = [ 10   1    1    1   0.12    0.12      10     0.6    0.6 ];

tvar  = {'2090s: \itp\rmCO_{2 AM}','2090s: RF','2090s: AF',...
         '2090s: \itp\rmCO_{2 AM} StDev','2090s: RF StDev','2090s: AF StDev',...
         '2090s - PI: \Delta\itp\rmCO_{2 AM}','2090s - PI: \DeltaRF','2090s - PI: \DeltaAF'}; 

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[1 1 6.6 4.25])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(3,3,[0 0],[0.05 .05],[0 0.03]);

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

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
    h = colorbar('east');    
    bn = crng1(n):cint(n):crng2(n);
    caxis([min(bn) max(bn)]);
    indc = [3 9];
    if sum(n == indc) >0
        cc = cmocean('-amp',length(bn));
    else
        cc = cmocean('amp',length(bn));
    end
    colormap(gca,cc);
    
    clear loc; loc = h.Position;    
    set(h,'position',[.01+loc(1) loc(2)  .5*loc(3) loc(4)]);
    set(h,'yaxislocation','right')
end


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for n = 1:length(let)
    axes(ha(n));hold on    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.1*xr,yylim(2)-0.125.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S3.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S3'])


%% Clean Up -----------------------------------

clear;close all;clc
