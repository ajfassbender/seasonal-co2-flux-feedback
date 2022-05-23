%% Figure 5
%
% DESCRIPTION:
%   Cumulative flux result maps for different chemical leverage configurations.
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


%%  Global Figures -----------------------------------

area  = getArea(-89.5:1:89.5,.5:359.5);
lat   = ncread([data_path '/Ensemble_Members/lat.nc'], 'lat');   
lon   = ncread([data_path '/Ensemble_Members/lon.nc'], 'lon'); 

fpath = [data_path '/Ensemble_Members'];

% Flux Differences
flux_l    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_l'); 
flux_nl   = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_nl');
flux_tl   = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_tl');
flux_bpl  = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_bpl');
d_l_nl    = flux_l   - flux_nl; 
d_tl_nl   = flux_tl  - flux_nl; 
d_bpl_nl  = flux_bpl - flux_nl;

var   = {'flux_l', 'd_l_nl', 'd_tl_nl', 'd_bpl_nl'};
crng1 = [-400 -100 -100 -100];
crng2 = [ 400  100  100  100];
cint  = [  80   20   20   20];
tvar  = {'Cumulative Flux (',...
         'Impact of Chem. Lev. (',...
         'Impact of T Lev. (',...
         'Impact of BP Lev. ('};

% Map projection limits
Slongs =[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats  = [0 90;-90 0;-90 0;0 90;-90 0;0 90];

x       = find(lon > 180);
llon    = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 12 6.5 3.8])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(2,2,[0 0],[0 .045],[0 0]);

% Average cumulative uptake across models
tflux_l     = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_l'); 
tflux_nl    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_nl');
tflux_tl    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_tl');
tflux_bpl   = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_bpl');
CU = [tflux_l; tflux_l - tflux_nl; tflux_l - tflux_bpl; tflux_l - tflux_tl];

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
    t1 = [char(tvar(n))  num2str(round(CU(n))) ' PgC)'];
    T= title({t1,'(mol C m^{-2})'},'fontweight','normal');set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h = colorbar;

    bn = crng1(n):cint(n):crng2(n);caxis([min(bn) max(bn)]);
    cc = cmocean('balance',length(bn));colormap(gca,cc);

    if n == 1
        set(h,'position',[0.072 0.58 0.0175 0.13])
    elseif n == 2
        set(h,'position',[0.571 0.58 0.0175 0.13])
    elseif n == 3
        set(h,'position',[0.072 0.1 0.0175 0.13])
    elseif n == 4
        set(h,'position',[0.571 0.1 0.0175 0.13])
    end
end


%% Subplot letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for n = 1:4
    axes(ha(n));hold on    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.05*xr,yylim(2)+0.06.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig=gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_5.pdf'],'Resolution',300)
savefig([fig_path '/Figure_5'])


%%  Table -----------------------------------

% Average cumulative uptake across models

flux_l     = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_l'); 
flux_nl    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_nl');
flux_tl    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_tl');
flux_bpl   = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_bpl');
flux_model = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_model'); 
x1 = [flux_l;flux_nl;flux_bpl;flux_tl;flux_model] % Pg C

% StDev of cumulative uptake (PgC) across runs;

flux_l_sd     = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_l_sd'); 
flux_nl_sd    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_nl_sd');
flux_tl_sd    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_tl_sd');
flux_bpl_sd   = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_bpl_sd');
flux_model_sd = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'tot_uptake_model_sd'); 
x2 = [flux_l_sd;flux_nl_sd;flux_bpl_sd;flux_tl_sd;flux_model_sd] % Pg C

vars    = {'Chem Lev','No Lev','BP Lev','T Lev','Model Output'};
delta   = x1 - x1(2);
delta_p = 100.*(x1 - x1(2))./x1(2);

t0 = cell2table(vars(:),'VariableNames',{'Cumulative Global Flux'});
t1 = array2table([x1 x2 delta delta_p],'VariableNames',{'EM Pg C' 'EM Pg C (StDev)' 'Uptake Diff' '% Diff'});
t  = [t0 t1]

% Differences relative to no-leverage
%  8.1% more uptake in reconstruction w/ leverage
% 18.4% less uptake in reconstruction w/ bp leverage only
% 10.2% more uptake in reconstruction w/ t leverage only

% Difference between reconstruction and model output
x1(end) - x1(1);
100.*(x1(end)- x1(1))./x1(1)
%  1.1% less uptake in model output relative to reconstruction w/ leverage

% Error Propagation for Leverage vs. No Leverage
v = (x2(1)^2 + x2(2)^2)^.5 
% --> 1.7284 PgC; uncertainty in the PgC difference (37 ± 2 Pg C)

val = 100 .* (x1(1) - x1(2))/ x1(2)
(val) .* ((v./(x1(1) - x1(2)))^2 + (x2(2)./x1(2))^2 ) ^.5
% --> 0.3801 % uncertainty in the percentage difference (8.1 ± 0.4 %)


%% North Pacific Percentage -----------------------------------

% regional boundaries [lat_min lat_max lon_min lon_mat];
r_ind = [ 30 42 140 240];

x1 = find(r_ind(1) <= lat & lat <= r_ind(2));
x2 = find(r_ind(3) <= lon & lon <= r_ind(4));

flux_l     = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_l'); 
flux_nl    = ncread([fpath '/EnsembleMean_C_Uptake.nc'],'uptake_nl');
f   = squeeze(nansum(nansum(flux_l(x2,x1,:)  .* transpose(area(x1,x2)))) .* 12 ./ (10^15));% Pg C/mo
fnl = squeeze(nansum(nansum(flux_nl(x2,x1,:) .* transpose(area(x1,x2)))) .* 12 ./ (10^15));% Pg C/mo

(f-fnl)./f % 14 percent

%% Clean up -----------------------------------

clear;close all;clc
