%% Figure 4
%
% DESCRIPTION:
%   1950s to 2090s ensemble mean pCO2 amplitude changes.
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

%% Figure Amplitude Changes -----------------------------------

fpath = [data_path '/Ensemble_Members/Ensemble_Mean_AmplitudeChanges_1950_2090.nc'];

co2_50s = ncread(fpath,'Aco2_totl_50s');   
dco2    = ncread(fpath,'Aco2_totl_90s') - co2_50s;   
dco2_p  = 100 .* dco2 ./ co2_50s;
x = find(dco2_p>200);dco2_p(x) = 200;
x = find(dco2_p<0);dco2_p(x)   = 0;

dco2_bp    = ncread(fpath,'Aco2_bp_90s')    - ncread(fpath,'Aco2_bp_50s');
dco2_bp_nl = ncread(fpath,'Aco2_bp_nl_90s') - ncread(fpath,'Aco2_bp_nl_50s');
dco2_bp_p  = 100 .* ((dco2_bp-dco2_bp_nl)./dco2_bp);
x = find(dco2_bp_p>100);dco2_bp_p(x) = 100;
x = find(dco2_bp_p<0);dco2_bp_p(x)   = 0;

dco2_t    = ncread(fpath,'Aco2_t_90s')    - ncread(fpath,'Aco2_t_50s');
dco2_t_nl = ncread(fpath,'Aco2_t_nl_90s') - ncread(fpath,'Aco2_t_nl_50s');
dco2_t_p  = 100 .* (dco2_t-dco2_t_nl)./dco2_t;
x = find(dco2_t_p>100);dco2_t_p(x) = 100;
x = find(dco2_t_p<0);dco2_t_p(x)   = 0;

lat  = ncread([data_path '/Ensemble_Members/lat.nc'],  'lat');   
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],  'lon'); 

figure(1);clf
set(gcf,'Units','Inches');
set(gcf,'position',[2 12 6.5 4.5])
set(gcf,'Resize','off')
clear ha;ha = tight_subplot(3,2,[0 0],[0 .045],[0 0]);

%map projection limits
Slongs=[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats= [0 90;-90 0;-90 0;0 90;-90 0;0 90];

var   = {'dco2', 'dco2_p', 'dco2_bp', 'dco2_bp_p', 'dco2_t', 'dco2_t_p'};
crng1 = [-30 -20 0 0 0 0];
crng2 = [300 200 300 100 300 100];
cint  = [30 20 30 10 30 10];
tvar  = {'\DeltaA-\itp\rmCO_2 (\muatm)',...
         '\DeltaA-\itp\rmCO_2 (%)',...
         '\DeltaA-\itp\rmCO_{2 BP} (\muatm)',...
         '\DeltaA-\itp\rmCO_{2 BP} Due to Chem. Lev. (%)',...
         '\DeltaA-\itp\rmCO_{2 T} (\muatm)',...
         '\DeltaA-\itp\rmCO_{2 T} Due to Chem. Lev. (%)'};

x = find(lon>180);
llon = lon;
llon(x) = lon(x)-360;
LON2 = [transpose(llon) .5];

for n = 1:length(var)    
    clear xvar
    eval(['xvar = ' char(var(n)) ';']);
    xvar = transpose(xvar);
    clear X;X = [xvar xvar(:,1)];
    
    axes(ha(n));hold on
    
    for l=1:6
        m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
        m_grid('fontsize',6,'xticklabels',[],'xtick',[],...
            'ytick',[],'yticklabels',[],'linest','-','color',[.5 .5 .5]);
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
    h = colorbar;
    
    bn = crng1(n):cint(n):crng2(n);
    if n ==1
        caxis([min(bn) max(bn)]);
        cd = cmocean('balance',length(bn));
        cc = [ cd(length(bn)/2 -1,:)  ; cmocean('amp',length(bn)-2)];
    elseif n ==2
        caxis([min(bn) max(bn)]);
        cd = cmocean('balance',length(bn));
        cc = [ cd(length(bn)/2 -1,:) ; cmocean('amp',length(bn)-2)];
    else
        caxis([min(bn) max(bn)]);
        cc = cmocean('amp',length(bn));
    end
    colormap(gca,cc);
    
    if n == 1
        set(h,'position',[0.073 0.685 0.0175 0.1])
    elseif n == 2
        set(h,'position',[0.572 0.685 0.0175 0.1])
    elseif n == 3
        set(h,'position',[0.073 0.3675 0.0175 0.1])
    elseif n == 4
        set(h,'position',[0.572 0.3675 0.0175 0.1])
    elseif n == 5
        set(h,'position',[0.073 0.049 0.0175 0.1])
    elseif n == 6
        set(h,'position',[0.572 0.049 0.0175 0.1])
    end
end

%% Subplot letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for n = 1:6
    axes(ha(n));hold on
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.05.*xr,yylim(2),let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',11,'fontweight','normal','fontname','Arial')
end

%% Save -----------------------------------

fig=gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_4.pdf'],'Resolution',300)
savefig([fig_path '/Figure_4'])

%% Clean up -----------------------------------

clear;close all;clc
