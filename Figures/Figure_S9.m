%% Figure S9
%
% DESCRIPTION:
%   Across ensemble member stdev of cumulative uptake.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   tight_subplot: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w 
%   m_map toolbox: https://www.eoas.ubc.ca/~rich/map.html
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

lon  = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 
lat  = ncread([data_path '/Ensemble_Members/lat.nc'],'lat');

l_sd   = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'uptake_l_sd');
nl_sd  = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'uptake_nl_sd');
bpl_sd = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'uptake_bpl_sd');
tl_sd  = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'uptake_tl_sd');

tot_l_sd   = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'tot_uptake_l_sd');
tot_nl_sd  = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'tot_uptake_nl_sd');
tot_bpl_sd = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'tot_uptake_bpl_sd');
tot_tl_sd  = ncread([data_path '/Ensemble_Members/EnsembleMean_C_Uptake.nc'],'tot_uptake_tl_sd');

%% Make Figure -----------------------------------

% Map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

% Cumulative mol/m2 
dnl  = (l_sd.^2 + nl_sd .^2).^.5;
dbpl = (l_sd.^2 + bpl_sd.^2).^.5;
dtl  = (l_sd.^2 + tl_sd .^2).^.5;

% Cumulative Pg C
tot_dnl  = (tot_l_sd.^2 + tot_nl_sd .^2).^.5
tot_dbpl = (tot_l_sd.^2 + tot_bpl_sd.^2).^.5
tot_dtl  = (tot_l_sd.^2 + tot_tl_sd .^2).^.5

var   = {'l_sd', 'dnl', 'dbpl', 'dtl'};
crng1 = [  0  0  0  0 ];
crng2 = [ 20 20 20 20 ];
cint  = [  2  2  2  2 ];
tvar  = {'Cumulative Flux (',...
         'Impact of Chem. Lev. (',...
         'Impact of T Lev. (',...
         'Impact of BP Lev. ('};

x       = find(lon > 180);
llon    = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

StDev = [tot_l_sd tot_dnl tot_dbpl tot_dtl];

figure(1);clf
set(gcf,'Units','inches');
set(gcf,'position',[1 1 4.33 3])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(2,2,[0 0],[0 .045],[0 0]);

for n = 1:length(var)
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
    t1 = [char(tvar(n))  char(177) num2str(round(StDev(n))) ' PgC)'];
    T  = title([{t1,'(mol C m^{-2})'}],'fontweight','normal');set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h  = colorbar;
        
    bn = crng1(n):cint(n):crng2(n);
    caxis([min(bn) max(bn)]);
    cc = cmocean('amp', length(bn));
    colormap(gca,cc);
    
    if n == 1
        set(h,'position',[0.072 0.61 0.0175 0.1])
    elseif n == 2
        set(h,'position',[0.571 0.61 0.0175 0.1])
    elseif n == 3
        set(h,'position',[0.072 0.132 0.0175 0.1])
    elseif n == 4
        set(h,'position',[0.571 0.132 0.0175 0.1])
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
exportgraphics(fig,[fig_path '/Figure_S9.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S9'])


%% Clean Up -----------------------------------

clear;close all;clc
