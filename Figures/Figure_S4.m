%% Figure S4
%
% DESCRIPTION:
%   Model output vs. reconstructed pCO2.
%   Ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   subtightplot: https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
%   m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
%   getArea: Computes surface area (m2) of each 1° latitide by 1° longitude horizontal grid
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

fpath = [data_path '/Ensemble_Members/m1/'];
vinfo = ncinfo([fpath 'ta.nc'],'ta');
dim = vinfo.Size;  

yr   = str2num(datestr(time,'yyyy'));
i90s = find(2090 <= yr & yr < 2100);
i20s = find(2020 <= yr & yr < 2030);
i50s = find(1950 <= yr & yr < 1960);

co2_2090s  = nanmean(ncread([fpath 'co2.nc'],'co2',[1 1 i90s(1)],[Inf Inf length(i90s)]),3);
dco2_2090s = co2_2090s - nanmean(ncread([fpath 'co2_totl.nc'],'co2_totl',[1 1 i90s(1)],[Inf Inf length(i90s)]),3);
co2_2020s  = nanmean(ncread([fpath 'co2.nc'],'co2',[1 1 i20s(1)],[Inf Inf length(i20s)]),3);
dco2_2020s = co2_2020s - nanmean(ncread([fpath 'co2_totl.nc'],'co2_totl',[1 1 i20s(1)],[Inf Inf length(i20s)]),3);
co2_1950s  = nanmean(ncread([fpath 'co2.nc'],'co2',[1 1 i50s(1)],[Inf Inf length(i50s)]),3);
dco2_1950s = co2_1950s - nanmean(ncread([fpath 'co2_totl.nc'],'co2_totl',[1 1 i50s(1)],[Inf Inf length(i50s)]),3);


%% Global Bias in Reconstructed pCO2 -----------------------------------

area = transpose(getArea(-89.5:1:89.5,0.5:1:359.5));
marea = nanmean(nanmean(area));
[nanmean(nanmean(abs(dco2_2090s.*area)))./marea nanstd(nanstd(dco2_2090s.*area))./marea]
[nanmean(nanmean(abs(dco2_2020s.*area)))./marea nanstd(nanstd(dco2_2020s.*area))./marea]
[nanmean(nanmean(abs(dco2_1950s.*area)))./marea nanstd(nanstd(dco2_1950s.*area))./marea]

%% Make Figure Lower Panels -----------------------------------

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 12 6.5 4])
set(gcf,'Resize','off')

loc     = {'KEO','SO'};
lat_ind = [130  32]; lat(lat_ind)
lon_ind = [160  10]; lon(lon_ind)
iyr     = [1950 2020 2090];

let  = {'(g)','(h)','(i)','(j)','(k)','(l)'};
titl = {'1950s','2020s','2090s'};

for n = 1:6
    subtightplot(6,3,n+12,[.04 .06],[0.06 .05],[0.15 0.03]);
    if n < 4
        q = find(iyr(n) <= yr & yr < iyr(n)+10);
        plot(time(q), squeeze(ncread([fpath 'co2_totl.nc'],'co2_totl',[lon_ind(1),lat_ind(1) q(1)],[1 1 length(q)])),'color',rgb('orange'),'linewidth',1);hold on
        plot(time(q), squeeze(ncread([fpath 'co2.nc'],'co2',[lon_ind(1),lat_ind(1) q(1)],[1 1 length(q)])),'color',rgb('dodgerblue'),'linewidth',1)
        plot([time(q(1)) time(q(end))], [0 0],'-k','linewidth',.5)        
        datetick('x','YY');xlim([time(q(1)) time(q(end))]);
        set(gca,'xticklabel',[])
        if n == 1
            ylim([200 400])
        elseif n == 2
            ylim([300 500])
        else
            ylim([700 1000])
        end
        title(titl(n))
    else
        q = find(iyr(n-3) <= yr & yr < iyr(n-3)+10);
        plot(time(q), squeeze(ncread([fpath 'co2_totl.nc'],'co2_totl',[lon_ind(2),lat_ind(2) q(1)],[1 1 length(q)])),'color',rgb('orange'),'linewidth',1);hold on
        plot(time(q), squeeze(ncread([fpath 'co2.nc'],'co2',[lon_ind(2),lat_ind(2) q(1)],[1 1 length(q)])),'color',rgb('dodgerblue'),'linewidth',1)
        plot([time(q(1)) time(q(end))], [0 0],'-k','linewidth',.5)
        datetick('x','YY');xlim([time(q(1)) time(q(end))]);
        if n == 4
            ylim([250 450])
        elseif n == 5
            ylim([380 480])
        else
            ylim([700 1000])
        end
    end
    grid on
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xx; xx = diff(xxlim);
    clear yy; yy = diff(yylim);
    text(xxlim(1)+0.02*xx,yylim(2)-0.4.*yy,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')

    if n == 1
        ty = ylabel({'NP','\itp\rmCO_2'});
        Z = ty.Position;
        set(ty,'Position',[Z(1)-250 Z(2) Z(3)])
        text(time(25),360,'Calc.','color',rgb('orange'))
        text(time(78),360,'Output','color',rgb('dodgerblue'))
    elseif n == 4
        ty = ylabel({'SO','\itp\rmCO_2'});
        R = ty.Position;
        set(ty,'Position',[Z(1)-250 R(2) R(3)])
    end
end


%% Make Figure Upper Panels -----------------------------------

% Map projection limits
Slongs =[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats  = [0 90;-90 0;-90 0;0 90;-90 0;0 90];


var = {'co2_1950s','co2_2020s','co2_2090s','dco2_1950s','dco2_2020s','dco2_2090s'};
crng1 = [250 350 775 -10 -10 -10];
crng2 = [400 500 925  10  10  10];
cint  = (crng2 - crng1)./10;
tvar  = {'1950s \itp\rmCO_2','Model Output','2020s \itp\rmCO_2','Model Output',...
         '2090s \itp\rmCO_2','Model Output','1950s \itp\rmCO_2','Output - Calc.',...
         '2020s \itp\rmCO_2','Output - Calc.','2090s \itp\rmCO_2','Output - Calc.'};

x = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2 = [transpose(llon) .5];

tind = 1;
for n = 1:length(var)
    clear xvar
    eval(['xvar = ' char(var(n)) ';']);
    xvar = transpose(xvar); % mol/m2/mo
    clear X;X = [xvar xvar(:,1)];
    
    if n < 4
        subtightplot(6,3,[n,n+3],[0 0],[0.05 .05],[0 0.03]);
    elseif n>3
        subtightplot(6,3,[n+3,n+6],[0 0],[0.05 .05],[0 0.03]);
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
        [xx,yy]=m_ll2xy(lon_ind(1),lat_ind(1)-90);
        if n == 4
            plot(xx,yy,'sw','markersize',5,'linewidth',1,'color',rgb('white'))
        else
            plot(xx,yy,'sk','markersize',5,'linewidth',1,'color',rgb('white'))
        end
        [xx,yy]=m_ll2xy(lon_ind(2),lat_ind(2)-90);
        plot(xx,yy,'sk','markersize',5,'linewidth',1,'color',rgb('white'))
    end    
    set(gca,'xlimmode','auto','ylimmode','auto');
    T= title({char(tvar(tind)), char(tvar(tind+1))});
    tind = tind + 2;
    set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h = colorbar('east');
    
    bn = crng1(n):cint(n):crng2(n);
    if n<4
        caxis([min(bn) max(bn)]);
        cc = cmocean('amp',length(bn));
        colormap(gca,cc);
    else
        caxis([min(bn) max(bn)]);
        cc = cmocean('balance',length(bn));
        colormap(gca,cc);
    end
    clear loc; loc = h.Position;    
    set(h,'position',[.01+loc(1) loc(2)  .5*loc(3) loc(4)]);
    set(h,'yaxislocation','right')
end

lon(lon_ind) %   9.5 159.5
lat(lat_ind) % -58.5  39.5


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for n = 1:6    
      if n < 4
        subtightplot(6,3,[n,n+3],[0 0],[0.05 .05],[0 0.03]);
    elseif n>3
        subtightplot(6,3,[n+3,n+6],[0 0],[0.05 .05],[0 0.03]);
      end
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.1*xr,yylim(2),let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S4.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S4'])


%% Clean Up -----------------------------------

clear;close all;clc
