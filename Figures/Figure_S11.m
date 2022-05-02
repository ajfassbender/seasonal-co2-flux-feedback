%% Figure S11
%
% DESCRIPTION:
%   Max/min pCO2 anomaly from pCO2_AM for chem-Lev (2090s; ensemble mean), no-Lev (2090s, ensemble mean),
%   and PI periods.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   tight_subplot: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w 
%   m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
%
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;
fig_path  = cd;


%% Ensemble Member 2090s Max/Min Calculations  -----------------------------------

vinfo = ncinfo([data_path '/Ensemble_Members/m1/ta.nc'],'ta');
dim = vinfo.Size; 

time  = ncread([data_path '/Ensemble_Members/time.nc'],'time');
yr    = str2num(datestr(time,'yyyy'));
uyr   = unique(yr); 
yri   = find(2090 <=  yr &  yr < 2100);
uyri  = find(2090 <= uyr & uyr < 2100);

% Pre-allocate
co2_mx     = NaN(dim(1),dim(2),30);
co2_mn     = NaN(dim(1),dim(2),30);
co2_nl_mx  = NaN(dim(1),dim(2),30);
co2_nl_mn  = NaN(dim(1),dim(2),30);

st_0 = now;
for i = 1:30
    clear co2_am co2 co2_nl
    disp(['Start ensemble member ' num2str(i)])   
    st_1 = now; 

    fpath = [data_path '/Ensemble_Members/m' num2str(i)];  

    co2_am = ncread([fpath '/co2_am.nc'],  'co2_am',  [1 1 uyri(1)],[Inf Inf length(uyri)]);
    co2    = ncread([fpath '/co2_totl.nc'],'co2_totl',[1 1  yri(1)],[Inf Inf length(yri)]);
    co2_nl = ncread([fpath '/co2_nl.nc'],  'co2_nl',  [1 1  yri(1)],[Inf Inf length(yri)]);

    co2_mx(:,:,i)    = max(nanmean(reshape(co2,dim(1),dim(2),12,[]),4),[],3) - nanmean(co2_am,3);
    co2_mn(:,:,i)    = min(nanmean(reshape(co2,dim(1),dim(2),12,[]),4),[],3) - nanmean(co2_am,3);
    co2_nl_mx(:,:,i) = max(nanmean(reshape(co2_nl,dim(1),dim(2),12,[]),4),[],3) - nanmean(co2_am,3);
    co2_nl_mn(:,:,i) = min(nanmean(reshape(co2_nl,dim(1),dim(2),12,[]),4),[],3) - nanmean(co2_am,3);

    disp(['Complete ensemble member ' num2str(i)])
    minutes((now - st_1)*24*60) 
end
disp('Total time:')
minutes((now - st_0)*24*60)


%% Control Run Max/Min Calculations  -----------------------------------

co2_av_pi = ncread([data_path '/Control_Run/climatology.nc'],'co2_mean');
co2_am_pi = ncread([data_path '/Control_Run/co2_am.nc'],'co2_am');
co2_pi_mx = max(co2_av_pi,[],3) - co2_am_pi;
co2_pi_mn = min(co2_av_pi,[],3) - co2_am_pi;


%% Make Figure -----------------------------------

lon  = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 
lat  = ncread([data_path '/Ensemble_Members/lat.nc'],'lat');

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

x       = find(lon > 180);
llon    = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

var   = {'nanmean(co2_mx,3)', 'nanmean(co2_mn,3)', 'nanmean(co2_nl_mx,3)', ...
         'nanmean(co2_nl_mn,3)','co2_pi_mx', 'co2_pi_mn'};
tvar  = {'2090s Max \itp\rmCO_{2 Anom.}', '2090s Min \itp\rmCO_{2 Anom.}','2090s Max \itp\rmCO_{2 \0 Anom.}',...
         '2090s Min \itp\rmCO_{2 \0 Anom.}', 'PI Max \itp\rmCO_{2 Anom.}','PI Min \itp\rmCO_{2 Anom.}'};

figure(1);clf
set(gcf,'Units','inches');
set(gcf,'position',[1 1 4.33 4.5])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(3,2,[0 0],[0 .045],[0 0]);

for n = 1:length(var)
    clear xvar 
    eval(['xvar = transpose(' char(var(n)) ');']);
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
    T = title([{char(tvar(n)),'(\muatm)'}],'fontweight','normal');set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h = colorbar;
        
    if n == 1 || n == 3 || n == 5
        bn = 0:25:150;
        cc = cmocean('amp', length(bn));colormap(gca,cc);
        caxis([min(bn) max(bn)]);        
    else n == 2 || n == 4 || n == 6   
        bn = -150:25:0;    
        cc = cmocean('-amp', length(bn));colormap(gca,cc);
        caxis([min(bn) max(bn)]);
    end
    
    if n == 1
        set(h,'position',[0.07 0.715 0.0175 0.0775])
    elseif n == 2
        set(h,'position',[0.571 0.715 0.0175 0.0775])
    elseif n == 3
        set(h,'position',[0.072 0.395 0.0175 0.0775])
    elseif n == 4
        set(h,'position',[0.571 0.395 0.0175 0.0775])
    elseif n == 5
        set(h,'position',[0.072 0.075 0.0175 0.0775])
    elseif n == 6
        set(h,'position',[0.571 0.075 0.0175 0.0775])
    end
end


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for n = 1:length(let)
    axes(ha(n));hold on
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.05*xr,yylim(2)+0.06.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S11.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S11'])


%% Clean Up -----------------------------------

clear;close all;clc
