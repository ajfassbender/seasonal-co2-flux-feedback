%% Figure 2
%
% DESCRIPTION:
%   Example of pCO2 Seasonal Cycle Changes at One Location
%   Use ensemble member 1 only
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   subtightplot: https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
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

lat  = ncread([data_path '/Ensemble_Members/lat.nc'],  'lat');   
lon  = ncread([data_path '/Ensemble_Members/lon.nc'],  'lon');  
time = ncread([data_path '/Ensemble_Members/time.nc'],'time'); 

% Specify Location
lats  = 120;lat(lats)
lons =  190;lon (lons)
clear lon lat

co2_totl  = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_totl.nc'], 'co2_totl', [lons lats 1],[1 1 Inf]));
co2_nl    = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_nl.nc'],   'co2_nl',   [lons lats 1],[1 1 Inf]));
co2_t     = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_t.nc'],    'co2_t',    [lons lats 1],[1 1 Inf]));
co2_bp    = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_bp.nc'],   'co2_bp',   [lons lats 1],[1 1 Inf])); 
co2_bp_nl = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_bp_nl.nc'],'co2_bp_nl',[lons lats 1],[1 1 Inf]));
co2_t_nl  = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_t_nl.nc'], 'co2_t_nl', [lons lats 1],[1 1 Inf]));
co2_atm   = squeeze(ncread([data_path '/Ensemble_Members/m2/co2_atm.nc'],  'co2_atm',  [lons lats 1],[1 1 Inf]));
flux_nl   = squeeze(ncread([data_path '/Ensemble_Members/m2/flux_keff_nl.nc'], 'flux_keff_nl',[lons lats 1],[1 1 Inf]));
flux_totl = squeeze(ncread([data_path '/Ensemble_Members/m2/flux_keff.nc'],    'flux_keff',   [lons lats 1],[1 1 Inf]));

co2_avg = ncread([data_path '/Ensemble_Members/m2/co2_am.nc'],'co2_am',[lons lats 1],[1 1 Inf]);
co2_am = reshape(repmat(transpose(squeeze(co2_avg)),12,1),[],1);
clear co2_avg


%% Make Figure -----------------------------------

% Specify three decades to highlight
yr     = str2num(datestr(time,'yyyy'));
t_ints = [transpose(find(1950<= yr & yr <1960));
          transpose(find(2040<= yr & yr <2050));
          transpose(find(2090<= yr & yr <2100))];
clear yr

figure(1);clf
set(gcf,'units','inches')
set(gcf,'position',[0 3 7 6.9])
set(gcf,'Resize','off')

% Create decade subplots
t_title = {'1950s','2040s','2090s'};
for j = 1:3    
    CO2_totl   = co2_totl(t_ints(j,:))  - co2_am(t_ints(j,:));
    CO2_nl     = co2_nl(t_ints(j,:))    - co2_am(t_ints(j,:));
    CO2_tl     = co2_t(t_ints(j,:))     - co2_am(t_ints(j,:));
    CO2_bpl    = co2_bp(t_ints(j,:))    - co2_am(t_ints(j,:));
    CO2_t_nl   = co2_t_nl(t_ints(j,:))  - co2_am(t_ints(j,:));
    CO2_bp_nl  = co2_bp_nl(t_ints(j,:)) - co2_am(t_ints(j,:));
    Time       = time(t_ints(j,:));
    
    % Create subplot (a)
    if j == 1
        subtightplot(5,3,1:3,.08,.03,[.12 .05])
        col = [rgb('honeydew'); rgb('cornsilk'); rgb('aliceblue')];
        for q = 1:3
            x = [t_ints(q,1) t_ints(q,1) t_ints(q,end) t_ints(q,end) t_ints(q,1)];
            y = [200 1100 1100 200 200];
            patch(time(x),y,col(q,:));hold on
            box on
        end
        y = co2_totl;
        plot(time,y,'k','linewidth',1.5,'color',rgb('orange'));hold on
        y = co2_nl;
        plot(time,y,'r','linewidth',1.5,'color',rgb('blueviolet'));hold on
        y = co2_atm;
        plot(time,y,'k','linewidth',2,'color',rgb('black'));hold on
        y = co2_am;
        plot(time,y,'-k','linewidth',1,'color',rgb('red'));hold on
        set(gca,'ytick',250:250:1000)
        dt = 1950:10:2100;
        m = 1;
        for ii = 1:1:length(dt)
            dd(m) =  datenum(['01/01/' num2str(dt(ii))],'mm/dd/yyyy');
            m = m+1;
        end
        set(gca,'xtick',dd);
        datetick('x','keepticks');xlim([datenum('1/1/1950') datenum('1/15/2100')])
        ylabel({'\itp\rmCO_2','(\muatm)'})
        text(time(450),880,'\itp\rmCO_{2}','color',rgb('orange'))
        text(time(450),650,'\itp\rmCO_{2 \0}','color',rgb('blueviolet'))
        text(time(660),880,'\itp\rmCO_{2 AM}','color',rgb('red'))
        text(time(660),650,'\itp\rmCO_{2 Atm.}','color',rgb('black'))
        ylim([200 1100]);xlim([time(1) time(end)]);grid on
        set(gca, 'Layer', 'Top')
        
        yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
        clear xr; xr = diff(xxlim);
        clear yr; yr = diff(yylim);
        text(xxlim(1)-.125.*xr,yylim(2)- .1*yr,'(a)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    end
    
    % Create subplot (k)
    if j == 1
        subtightplot(5,3,13:15,.08,[.05 .03],[.12 .05])
        col = [rgb('honeydew'); rgb('cornsilk'); rgb('aliceblue')];
        for q = 1:3
            x = [t_ints(q,1) t_ints(q,1) t_ints(q,end) t_ints(q,end) t_ints(q,1)];
            y = [-25.1 25.1 25.1 -25.1 -25.1];
            patch(time(x),y,col(q,:));hold on
            box on
        end
        y = 1000.*flux_totl;
        plot(time,y,'k','linewidth',1.5,'color',rgb('orange'));hold on
        y = 1000.*flux_nl;
        plot(time,y,'r','linewidth',1.5,'color',rgb('blueviolet'));hold on
        plot([time(1) time(end)],[0 0],'k','linewidth',1);hold on
        set(gca,'ytick',-24:12:24);
        dt = 1950:10:2100;
        m = 1;
        for ii = 1:1:length(dt)
            dd(m) =  datenum(['01/01/' num2str(dt(ii))],'mm/dd/yyyy');
            m = m+1;
        end
        set(gca,'xtick',dd);
        datetick('x','keepticks');xlim([datenum('1/1/1950') datenum('1/15/2100')])
        ylabel({'CO_2 Flux','(mmol m^{-2} mo^{-1})'})
        ylim([-24 24]);xlim([time(1) time(end)]);grid on
        set(gca, 'Layer', 'Top')        
        yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
        clear xr; xr = diff(xxlim);
        clear yr; yr = diff(yylim);
        text(xxlim(1)-.125.*xr,yylim(2)+ .1*yr,'(k)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    end
    
    % Create thermal pCO2 component subplots
    h =  subtightplot(5,3,j+3,.03,.05,[.12 .05]);
    plot(Time,CO2_tl,'k','linewidth',1.5,'color',rgb('orange'));hold on
    plot(Time,CO2_t_nl,'-r','linewidth',1.5,'color',rgb('blueviolet'));hold on
    plot([Time(1) Time(end)],[0 0],'-k')
    set(gca,'xtick',Time(1):2.*365:Time(end))
    datetick('x','yy','keepticks');
    xlim([Time(1) Time(end)])
    title([char(t_title(j))])
    ylim([-275 275])
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    if j+3 ==4
        ylabel({'\itp\rmCO_{2 T} Anom.','(\muatm)'})
        set(h,'Color',col(1,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(b)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    elseif j+3 == 5
        set(gca,'yticklabel',[])
        set(h,'Color',col(2,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(c)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    elseif j+3 == 6
        set(gca,'yticklabel',[])
        set(h,'Color',col(3,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(d)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    end
    box on; grid on
    set(gca,'xticklabel',[])
    
    % Create biphysical pCO2 component subplots
    h =  subtightplot(5,3,j+6,.03,.05,[.12 .05]);
    plot(Time,CO2_bpl,'k','linewidth',1.5,'color',rgb('orange'));hold on
    plot(Time,CO2_bp_nl,'-r','linewidth',1.5,'color',rgb('blueviolet'));hold on
    plot([Time(1) Time(end)],[0 0],'-k')
    set(gca,'xtick',Time(1):2.*365:Time(end))
    datetick('x','yy','keepticks');
    xlim([Time(1) Time(end)])
    ylim([-275 275])
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    if  j+6== 7
        ylabel({'\itp\rmCO_{2 BP} Anom.','(\muatm)'})
        set(h,'Color',col(1,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(e)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    elseif j+6 == 8
        set(gca,'yticklabel',[])
        set(h,'Color',col(2,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(f)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    elseif j+6 == 9
        set(gca,'yticklabel',[])
        set(h,'Color',col(3,:));
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(g)','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8)
    end
    box on; grid on
    set(gca,'xticklabel',[])
    
    % Create total pCO2 component subplots
    h =  subtightplot(5,3,j+9,.03,.05,[.12 .05]);
    plot(Time,CO2_totl,'k','linewidth',1.5,'color',rgb('orange'));hold on
    plot(Time,CO2_nl,'-r','linewidth',1.5,'color',rgb('blueviolet'));hold on
    plot([Time(1) Time(end)],[0 0],'-k')
    set(gca,'xtick',Time(1):2.*365:Time(end))
    datetick('x','yy','keepticks');
    xlim([Time(1) Time(end)])
    ylim([-275 275])
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);
    if j+9==10
        ylabel({'\itp\rmCO_{2} Anom.','(\muatm)'})
        set(h,'Color',col(1,:));        
        text(xxlim(1)-.1*xr,yylim(2)+ .01*yr,'(h)','VerticalAlignment','bottom',...
            'HorizontalAlignment','left','fontsize',8)
    elseif j+9 == 11
        set(gca,'yticklabel',[])
        set(h,'Color',col(2,:));       
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(i)','VerticalAlignment','bottom',...
            'HorizontalAlignment','left','fontsize',8)
    elseif j+9 == 12
        set(gca,'yticklabel',[])
        set(h,'Color',col(3,:));        
        text(xxlim(1)-.1.*xr,yylim(2)+ .01*yr,'(j)','VerticalAlignment','bottom',...
            'HorizontalAlignment','left','fontsize',8)
    end
    box on; grid on
end

%% Save -----------------------------------

fig = gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_2.pdf'],'Resolution',300)
savefig([fig_path '/Figure_2'])

%% Clean up -----------------------------------

clear;close all;clc
