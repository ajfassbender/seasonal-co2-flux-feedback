%% Figure S8
%
% DESCRIPTION:
%   Ensemble mean 1950s vs. 2090s seasonal cycle amplitudes/
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   subtightplot: https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
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

fpath = [data_path '/Ensemble_Members/'];
vars = {'Adic_50s','Adic_90s','Asst_50s','Asst_90s','Awind_50s','Awind_90s','Ata_50s','Ata_90s','Aco2_totl_50s','Aco2_totl_90s','mld_max_50s','mld_max_90s'};
for i = 1:length(vars)
    uvar = char(vars(i));
    x = ncread([fpath 'Ensemble_Mean_AmplitudeChanges_1950_2090.nc'], uvar);
    eval([uvar ' = x;']);
end


%% Make Figure -----------------------------------

lon = ncread([data_path '/Ensemble_Members/lon.nc'],'lon'); 
lat = ncread([data_path '/Ensemble_Members/lat.nc'],'lat');

% Map projection limits
Slongs =[-100.5 46.5;-76.5 25.5; 25.5 142.5;46.5 97.5;142.5 295;97.5 295];
Slats  = [0 90;-90 0;-90 0;0 90;-90 0;0 90];

ddic  = Adic_90s  - Adic_50s;
dsst  = Asst_90s  - Asst_50s;
dwind = Awind_90s - Awind_50s;
dta   = Ata_90s   - Ata_50s;
dmld  = mld_max_90s - mld_max_50s;
dco2  = Aco2_totl_90s - Aco2_totl_50s;

var   = { 'Adic_50s','Adic_90s','ddic',...
          'Ata_50s', 'Ata_90s', 'dta',...
          'Asst_50s','Asst_90s','dsst',...
          'Awind_50s','Awind_90s','dwind',...
          'Aco2_totl_50s','Aco2_totl_90s','dco2'...
          'mld_max_50s','mld_max_90s','dmld'};

crng1 = [  0     0   -40     0     0   -15      0     0   -2      0    0   -1.25     0     0   -25    0    0  -100];
crng2 = [125   125    40    75    75    15     15    15    2     10   10    1.25   200   200   125  250  250    50];
cint  = [12.5   12.5   0.8   7.5   7.5   0.15   1.5   1.5  0.4    1    1    0.3     20    20    15   25   25    15];

tvar  = {'1950s A-DIC','2090s A-DIC','\DeltaA-DIC',...
         '1950s A-TA','2090s A-TA','\DeltaA-TA',...
         '1950s A-SST','2090s A-SST','\DeltaA-SST',...
         '1950s A-Wind Speed','2090s A-Wind Speed','\DeltaA-Wind Speed',...
         '1950s A-\itp\rmCO_2','2090s A-\itp\rmCO_2','\DeltaA-\itp\rmCO_2',...
         '1950s Max MLD','2090s Max MLD','\DeltaMax MLD'};

unit = {'(\mumol kg^{-1})','(\mumol kg^{-1})','(\mumol kg^{-1})',...
        '(\mumol kg^{-1})','(\mumol kg^{-1})','(\mumol kg^{-1})',...
        '(\circC)','(\circC)','(\circC)',...
        '(m s^{-1})','(m s^{-1})','(m s^{-1})',...
        '(\muatm)','(\muatm)','(\muatm)',...
        '(m)','(m)','(m)'};


x       = find(lon > 180);
llon    = lon;
llon(x) = lon(x)-360;
LON2    = [transpose(llon) .5];

figure(1);clf
set(gcf,'Units','inches');
set(gcf,'position',[2 12 6.5 7])
set(gcf,'Resize','off')

for n = 1:length(var) 
    subtightplot(6,3,n,[0 0],[0.05 0.05],[0 0.05]);
    
    % grab data
    clear xvar 
    eval(['xvar = ' char(var(n)) ';']);
    xvar = transpose(xvar);
    clear X;X = [xvar xvar(:,1)];
        
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
    T = title({char(tvar(n)), char(unit(n))},'fontweight','normal');set(T,'position',[1.0000 1.0098+.08 -5.0000e+15])
    h = colorbar;
        
    bn = crng1(n):cint(n):crng2(n);
    caxis([min(bn) max(bn)]);
    if n == 3 || n == 6 || n == 9 || n == 12 || n == 15 || n == 18
        cc = cmocean('balance', length(bn),'pivot',0);
    else
        cc = cmocean('amp', length(bn));
    end
    colormap(gca,cc);

    clear loc; loc = h.Position;    
    set(h,'position',[.045+loc(1) loc(2)  .7*loc(3) loc(4)]);
    set(h,'yaxislocation','right')
end


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)',...
       '(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)'};
for n = 1:length(let)
    subtightplot(6,3,n,[0 0],[0.05 0.05],[0 0.05]);hold on    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    text(xxlim(1)+0.1*xr,yylim(2),let(n),'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',10,'fontweight','normal','fontname','Arial')
end


%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_S8.pdf'],'Resolution',300)
savefig([fig_path '/Figure_S8'])


%% Clean Up -----------------------------------

clear;close all;clc
