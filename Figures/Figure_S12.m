%% Figure S12
%
% DESCRIPTION:
%   2090s max/min delta pCO2 Sea-Ar for chem-leverage and no-leverage conditions.
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


%% Load Data -----------------------------------

fpath = [data_path '/Ensemble_Members/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'];
dco2    = ncread([fpath], 'delta_co2');
dco2_nl = ncread([fpath], 'delta_co2_nl');

% average across ensemble
mx    = transpose(nanmean(max(dco2,[],3),3));
mn    = transpose(nanmean(min(dco2,[],3),3));
mx_nl = transpose(nanmean(max(dco2_nl,[],3),3));
mn_nl = transpose(nanmean(min(dco2_nl,[],3),3));


%% Make Figure -----------------------------------

lon  = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 
lat  = ncread([data_path '/Ensemble_Members/lat.nc'],'lat');

% Map projection limits
Slongs =[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats  = [0 90;-90 0;-90 0;0 90;-90 0;0 90];

var   = {'mx', 'mn', 'mx_nl', 'mn_nl'};
crng1 = [ -25   -150   -25 -150 ];
crng2 = [ 150     25   150   25 ];
cint  = [  25     25    25   25 ];
tvar  = {'2090s Max \Delta\itp\rmCO_{2 Sea-Air}',...
         '2090s Min \Delta\itp\rmCO_{2 Sea-Air}',...
         '2090s Max \Delta\itp\rmCO_{2 \0 Sea-Air}',...
         '2090s Min \Delta\itp\rmCO_{2 \0 Sea-Air}'};

x       = find(lon > 180);
llon    = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

figure(1);clf
set(gcf,'Units','inches');
set(gcf,'position',[1 1 4.33 3])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(2,2,[0 0],[0 .045],[0 0]);

for n = 1:length(var)
    clear xvar 
    eval(['xvar = ' char(var(n)) ';']);
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
        
    bn = crng1(n):cint(n):crng2(n);
    caxis([min(bn) max(bn)]);
    if n == 1 || n == 3
    cc = cmocean('amp', length(bn));colormap(gca,cc);
    else    	
    cc = cmocean('-amp', length(bn));colormap(gca,cc);
    end
    
    if n == 1
        set(h,'position',[0.072 0.598 0.0175 0.11])
    elseif n == 2
        set(h,'position',[0.571 0.598 0.0175 0.11])
    elseif n == 3
        set(h,'position',[0.072 0.12 0.0175 0.11])
    elseif n == 4
        set(h,'position',[0.571 0.12 0.0175 0.11])
    end
end


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)'};
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
exportgraphics(fig,[fig_path '/Figure_S12.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S12'])


%% Clean Up -----------------------------------

clear;close all;clc
