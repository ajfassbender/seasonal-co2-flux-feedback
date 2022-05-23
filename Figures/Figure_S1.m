%% Figure S1
%
% DESCRIPTION:
%   Carbonate system internal consistency figure.
%   Use ensemble member 1 only
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


%% Compute 2090s pCO2 and DIC Differences -----------------------------------

time = ncread([data_path '/Ensemble_Members/time.nc'],'time'); 
yr   = str2num(datestr(time,'yyyy'));
ind  = find(2090 <= yr & yr < 2100);

fpath = [data_path '/Ensemble_Members/m1/'];
co2 = ncread([fpath 'co2.nc'], 'co2',[1 1 ind(1)],[Inf Inf length(ind)]);% model output
dic = ncread([fpath 'dic.nc'], 'dic',[1 1 ind(1)],[Inf Inf length(ind)]);% model output
co2_calc = ncread([fpath 'co2_calc.nc'], 'co2_calc',[1 1 ind(1)],[Inf Inf length(ind)]);% calculated
dic_calc = ncread([fpath 'dic_calc.nc'], 'dic_calc',[1 1 ind(1)],[Inf Inf length(ind)]);% calculated
dco2 = nanmean(co2,3) - nanmean(co2_calc,3);
ddic = nanmean(dic,3) - nanmean(dic_calc,3);


%% Make Figure -----------------------------------

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 12 6.5 1.75])
set(gcf,'Resize','off')
%clear ha;ha = tight_subplot(1,3,[0 0],[0.05 .05],[0 0]);

lat  = ncread([data_path '/Ensemble_Members/lat.nc'],  'lat');   
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],  'lon');

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

% Set figure titles and colorbar ranges
var   = {'dco2','ddic'};
crng1 = [-40 -10];
crng2 = [ 40  10];
cint  = [  8   2];
tvar  = {'2090s: \itp\rmCO_{2} Model - Calc.','2090s: DIC Model - Calc.'};

x    = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

% Set lat/lon for BATS
llat =  31 + 40/60;
llon = 360 - 64 + (10/60);
[~,bats_lon] = min(abs(lon - llon));
[~,bats_lat] = min(abs(lat - llat));

% Figure Panels 1 and 2
for n = 1:length(var)    
    eval(['xvar = ' char(var(n)) ';']);
    xvar = transpose(xvar);
    X    = [xvar xvar(:,1)];
    
   if n ==1 
    subtightplot(2,3,[1,4],[0 0],[0.05 .05],[0 0]);
   else
    subtightplot(2,3,[2,5],[0 0],[0.05 .05],[0 0]);
   end
    
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
            [xx,yy]=m_ll2xy(lon(bats_lon)-360,lat(bats_lat));
            plot(xx,yy,'sw','markersize',5,'linewidth',1)
        end
    
    set(gca,'xlimmode','auto','ylimmode','auto');
    T = title(tvar(n));ty = T.Position;
    if n == 2
        set(T,'position',[ty(1) ty(2)+0.2*ty(2) ty(3)])
    end
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

% Figure Panel 3

% Load Data from BATS Location
co2 = squeeze(ncread([fpath 'co2.nc'], 'co2',[bats_lon bats_lat 1],[1 1 Inf]));
dic = squeeze(ncread([fpath 'dic.nc'], 'dic',[bats_lon bats_lat 1],[1 1 Inf]));
co2_calc = squeeze(ncread([fpath 'co2_calc.nc'], 'co2_calc',[bats_lon bats_lat 1],[1 1 Inf]));
dic_calc = squeeze(ncread([fpath 'dic_calc.nc'], 'dic_calc',[bats_lon bats_lat 1],[1 1 Inf]));

subtightplot(2,3,3,[.2 .1],[0.1 .15],[0.03 .03]);
plot(time,co2-co2_calc,'r','color',cc(10,:));hold on
title('\itp\rmCO_{2} Model - Calc.')
set(gca,'xticklabel',[]);ylabel('\muatm')
set(gca,'ytick',0:25:50,'xtick',[datenum('1/1/2000') datenum('1/1/2050') datenum('1/1/2100')])
xlim([time(1) time(end)]);grid on

subtightplot(2,3,6,[.2 .1],[0.15 .15],[0.03 .03]);
plot(time,dic-dic_calc,'b','color',cc(3,:));hold on
title('DIC Model - Calc.');ylabel('\mumol kg^{-1}')
set(gca,'ytick',-10:5:0);ylim([-10 0])
datetick('x');xlim([time(1) time(end)]);grid on


%% Subplot Letters -----------------------------------

subtightplot(2,3,[1,4],[0 0],[0.05 .05],[0 0]);
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)+0.1*xr,yylim(2)-0.125.*yr,'(a)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')

subtightplot(2,3,[2,5],[0 0],[0.05 .05],[0 0]);
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)+0.1*xr,yylim(2)-0.125.*yr,'(b)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')

subtightplot(2,3,3,[.2 .1],[0.1 .15],[0.03 .03]);
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)-0.2*xr,yylim(2)+0.2.*yr,'(c)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')

subtightplot(2,3,6,[.2 .1],[0.15 .15],[0.03 .03]);
yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
clear xr; xr = diff(xxlim);
clear yr; yr = diff(yylim);
text(xxlim(1)-0.2*xr,yylim(2)+0.2.*yr,'(d)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S1.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S1'])


%% Clean Up -----------------------------------

clear;close all;clc
