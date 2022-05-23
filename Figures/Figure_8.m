%% Figure 8
%
% DESCRIPTION:
%   Regional results: 2090s seasonal cycles and 2090s and PI sensitivity terms.
%
% USER INPUT:
%   data_path: directory containing reprocessed .nc files
%   fig_path:  directory where figure will be saved
%
% FUNCTIONS CALLED:
%   subtightplot: https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
%   area_weighted_mean: computes area weighted mean
%   rgb: https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2/
% 
% AUTHOR:
%   A. J. Fassbender (NOAA-PMEL): andrea.j.fassbender@noaa.gov
%
% DATE: MAY 2, 2022


%% Set Directory Paths -----------------------------------

disp('set your paths')
data_path = cd;
fig_path  = cd;

%% Set Regional Bounaries -----------------------------------

% REGION = [Lon_min Lon_max Lat_min Lat_max]

NASP = [315 350 140 155];
NPSP = [150 185 132 145];
NPST = [145 225 120 130];
Eq   = [210 275  80 100];
SO   = [  1 360  40  60];

reg = [NASP; NPSP; NPST; Eq; SO];
leg = {'NA Subpolar','NP Subpolar','NP Subtropic','Eq. Pacific','S. Ocean'};

%% Area Weight Ensemble Mean Seasonal Cycles Over Regions -----------------------------------

fpath = [data_path '/Ensemble_Members'];
lat   = ncread([fpath '/lat.nc'],  'lat');   
lon   = ncread([fpath '/lon.nc'],  'lon'); 

vars  = {'delta_co2','delta_co2_nl','wind_speed','co2_flux','co2_flux_nl','delta_co2_bp','delta_co2_t',...
         'delta_co2_bp_nl','delta_co2_t_nl','rf','af','co2_am','rf_pi','af_pi','co2_am_pi'};

for i = 1:length(leg)
  clear lat_1 lon_1
  lon_1(:,1) = reg(i,1) : reg(i,2);
  lat_1(:,1) = reg(i,3) : reg(i,4);
  for j = 1:length(vars)
    clear x xx x_out
    uvar = char(vars(j));
    if j>9
      x   = ncread([fpath '/Figure8_2090s_Sensitivity_Terms.nc'],uvar);
      xx  = x(lon_1,lat_1,:);
      x_out = area_weighted_mean(xx,lat(lat_1),lon_1);
      eval([uvar '_awm(:,i) = x_out;'])
      else
        for mm = 1:12
          x   = ncread([fpath '/Figure8_2090s_Ensemble_Mean_SeasCycles.nc'],uvar,[1 1 mm], [Inf Inf 1]);  
          xx  = x(lon_1,lat_1,:);
          x_out(mm) = area_weighted_mean(xx,lat(lat_1),lon_1);  
        end
        eval([uvar '_awm(:,i) = x_out;'])
      end
  end
end
clear lon_1 lat_1 lon_2 lat_2 x xx vars i j uvar x_out


%% Make Figure -----------------------------------

figure(1);clf
set(gcf,'units','inches')
set(gcf,'position',[1 1 7 5.5])
set(gcf,'Resize','off')

%reg = [NPSP; PST; NASP; AST; Eq;  IO; SO];

x = transpose(1:12);

% NASP -----------------------------------
n = 1;
subtightplot(5,3,1,[0.03 0.11],[0.05 .05],[0.11 0.08]);
plot(x,delta_co2_t_awm(:,n),'linewidth',4,'color',rgb('Orange'));hold on
plot(x,delta_co2_t_nl_awm(:,n),'linewidth',2,'color',rgb('OrangeRed'));
plot(x,delta_co2_bp_awm(:,n),'linewidth',4,'color',rgb('LimeGreen'));
plot(x,delta_co2_bp_nl_awm(:,n),'linewidth',2,'color',rgb('ForestGreen'));
ylabel({char(leg(n)),'(\muatm)'})
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
title({'\Delta\itp\rmCO_2'})
ylim([-275 200])
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:100:300)

subtightplot(5,3,2,[0.03 0.11],[0.05 .05],[0.11 0.08]);
yyaxis left
plot(x,delta_co2_awm(:,n),'-k','linewidth',4,'color',rgb('LightSkyBlue'));hold on
plot(x,delta_co2_nl_awm(:,n),'-k','linewidth',2,'color',rgb('DodgerBlue'));
plot(x,co2_flux_awm(:,n).*100,'-k','linewidth',4,'color',rgb('silver'));hold on
plot(x,co2_flux_nl_awm(:,n).*100,'-k','linewidth',2,'color',rgb('black'));
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
T2 = title({'\Delta\itp\rmCO_2 & F\times100 | Wind Speed'});
ylabel({'\muatm &','mol m^{-2} mo^{-1}'},'color','k')
ylim([-150 25])
ax = gca;
ax.YColor = rgb('black');
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:50:300)

yyaxis right 
plot(x,wind_speed_awm(:,n),'-k','linewidth',3,'color',rgb('deeppink'));hold on
ylim([5 15])
ax = gca;
ax.YColor = rgb('deeppink');
ylabel('m s^{-1}')


% NPST -----------------------------------
n = 2;
subtightplot(5,3,4,[0.03 0.11],[0.05 .05],[0.11 0.08]);
plot(x,delta_co2_t_awm(:,n),'linewidth',4,'color',rgb('Orange'));hold on
plot(x,delta_co2_t_nl_awm(:,n),'linewidth',2,'color',rgb('OrangeRed'));
plot(x,delta_co2_bp_awm(:,n),'linewidth',4,'color',rgb('LimeGreen'));
plot(x,delta_co2_bp_nl_awm(:,n),'linewidth',2,'color',rgb('ForestGreen'));
ylabel({char(leg(n)),'(\muatm)'})
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylim([-275 200])
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:100:300)

subtightplot(5,3,5,[0.03 0.11],[0.05 .05],[0.11 0.08]);
yyaxis left
plot(x,delta_co2_awm(:,n),'-k','linewidth',4,'color',rgb('LightSkyBlue'));hold on
plot(x,delta_co2_nl_awm(:,n),'-k','linewidth',2,'color',rgb('DodgerBlue'));
plot(x,co2_flux_awm(:,n).*100,'-k','linewidth',4,'color',rgb('silver'));hold on
plot(x,co2_flux_nl_awm(:,n).*100,'-k','linewidth',2,'color',rgb('black'));
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylabel({'\muatm &','mol m^{-2} mo^{-1}'},'color','k')
ylim([-150 100])
ax = gca;
ax.YColor = rgb('black');
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:50:300)

yyaxis right 
plot(x,wind_speed_awm(:,n),'-k','linewidth',3,'color',rgb('deeppink'));hold on
ylim([5 15])
ax = gca;
ax.YColor = rgb('deeppink');
ylabel('m s^{-1}')


% NAST -----------------------------------
n = 3;
subtightplot(5,3,7,[0.03 0.11],[0.05 .05],[0.11 0.08]);
plot(x,delta_co2_t_awm(:,n),'linewidth',4,'color',rgb('Orange'));hold on
plot(x,delta_co2_t_nl_awm(:,n),'linewidth',2,'color',rgb('OrangeRed'));
plot(x,delta_co2_bp_awm(:,n),'linewidth',4,'color',rgb('LimeGreen'));
plot(x,delta_co2_bp_nl_awm(:,n),'linewidth',2,'color',rgb('ForestGreen'));
YS1 = ylabel({char(leg(n)),'(\muatm)'});
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylim([-275 200])
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:100:200)

subtightplot(5,3,8,[0.03 0.11],[0.05 .05],[0.11 0.08]);
yyaxis left
plot(x,delta_co2_awm(:,n),'-k','linewidth',4,'color',rgb('LightSkyBlue'));hold on
plot(x,delta_co2_nl_awm(:,n),'-k','linewidth',2,'color',rgb('DodgerBlue'));
plot(x,co2_flux_awm(:,n).*100,'-k','linewidth',4,'color',rgb('silver'));hold on
plot(x,co2_flux_nl_awm(:,n).*100,'-k','linewidth',2,'color',rgb('black'));
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
YS2 = ylabel({'\muatm &','mol m^{-2} mo^{-1}'},'color','k')
ylim([-100 125])
ax = gca;
ax.YColor = rgb('black');
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-300:50:200)

yyaxis right 
plot(x,wind_speed_awm(:,n),'-k','linewidth',3,'color',rgb('deeppink'));hold on
ylim([5 15])
ax = gca;
ax.YColor = rgb('deeppink');
ylabel('m s^{-1}')


% Eq -----------------------------------
n = 4;
subtightplot(5,3,10,[0.03 0.11],[0.05 .05],[0.11 0.08]);
plot(x,delta_co2_t_awm(:,n),'linewidth',4,'color',rgb('Orange'));hold on
plot(x,delta_co2_t_nl_awm(:,n),'linewidth',2,'color',rgb('OrangeRed'));
plot(x,delta_co2_bp_awm(:,n),'linewidth',4,'color',rgb('LimeGreen'));
plot(x,delta_co2_bp_nl_awm(:,n),'linewidth',2,'color',rgb('ForestGreen'));
ys1 = ylabel({char(leg(n)),'(\muatm)'});
set(ys1,'position',[YS1.Position(1) ys1.Position(2:3)])
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylim([-50 50])
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-50:25:50)

subtightplot(5,3,11,[0.03 0.11],[0.05 .05],[0.11 0.08]);
yyaxis left
plot(x,delta_co2_awm(:,n),'-k','linewidth',4,'color',rgb('LightSkyBlue'));hold on
plot(x,delta_co2_nl_awm(:,n),'-k','linewidth',2,'color',rgb('DodgerBlue'));
plot(x,co2_flux_awm(:,n).*100,'-k','linewidth',4,'color',rgb('silver'));hold on
plot(x,co2_flux_nl_awm(:,n).*100,'-k','linewidth',2,'color',rgb('black'));
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ys2 = ylabel({'\muatm &','mol m^{-2} mo^{-1}'},'color','k');
set(ys2,'position',[YS2.Position(1) ys2.Position(2:3)])
ylim([-15 15])
ax = gca;
ax.YColor = rgb('black');
set(gca,'xtick',1:1:12,'xticklabel',[],'ytick',-50:5:50)

yyaxis right 
plot(x,wind_speed_awm(:,n),'-k','linewidth',3,'color',rgb('deeppink'));hold on
ylim([5 15])
ax = gca;
ax.YColor = rgb('deeppink');
ylabel('m s^{-1}')


% Southern Ocean -----------------------------------
n = 5;
subtightplot(5,3,13,[0.03 0.11],[0.05 .05],[0.11 0.08]);
plot(x,delta_co2_t_awm(:,n),'linewidth',4,'color',rgb('Orange'));hold on
plot(x,delta_co2_t_nl_awm(:,n),'linewidth',2,'color',rgb('OrangeRed'));
plot(x,delta_co2_bp_awm(:,n),'linewidth',4,'color',rgb('LimeGreen'));
plot(x,delta_co2_bp_nl_awm(:,n),'linewidth',2,'color',rgb('ForestGreen'));
ylabel({char(leg(n)),'(\muatm)'})
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylim([-150 150])
set(gca,'xtick',1:1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'ytick',-150:75:150)

subtightplot(5,3,14,[0.03 0.11],[0.05 .05],[0.11 0.08]);
yyaxis left
plot(x,delta_co2_awm(:,n),'-k','linewidth',4,'color',rgb('LightSkyBlue'));hold on
plot(x,delta_co2_nl_awm(:,n),'-k','linewidth',2,'color',rgb('DodgerBlue'));
plot(x,co2_flux_awm(:,n).*100,'-k','linewidth',4,'color',rgb('silver'));hold on
plot(x,co2_flux_nl_awm(:,n).*100,'-k','linewidth',2,'color',rgb('black'));
plot([1 12],[0 0],'-k')
xlim([1 12]);grid on
ylabel({'\muatm &','mol m^{-2} mo^{-1}'},'color','k')
ylim([-75 25])
ax = gca;
ax.YColor = rgb('black');
set(gca,'xtick',1:1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'ytick',-150:25:150)

yyaxis right 
plot(x,wind_speed_awm(:,n),'-k','linewidth',3,'color',rgb('deeppink'));hold on
ylim([5 15])
ax = gca;
ax.YColor = rgb('deeppink');
ylabel('m s^{-1}')


%% Sensitivity Term Changes -----------------------------------

for i = 1:5
  subtightplot(5,3,i*3,[0.03 0.11],[0.05 .05],[0.11 0.08]);
  yyaxis left

  H = bar(2,rf_awm(i),'FaceColor','flat','BarWidth',.8,'EdgeColor','none');hold on
  H.CData = rgb('LimeGreen');
  H = bar(1,rf_pi_awm(i),'FaceColor','flat','BarWidth',.8,'EdgeColor','none');hold on
  H.CData = rgb('ForestGreen');

  H = bar(5,af_awm(i),'FaceColor','flat','BarWidth',.8,'EdgeColor',('none'));hold on
  H.CData = rgb('LimeGreen');
  H = bar(4,af_pi_awm(i),'FaceColor','flat','BarWidth',.8,'EdgeColor',('none'));hold on
  H.CData = rgb('ForestGreen');

  grid on
  ylim([-20 20])
  set(gca,'ytick',-20:10:20)
  ax = gca;
  ax.YColor = rgb('ForestGreen'); 
  set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',[])

  if i == 1 
    text(0.6,-5,'PI')
    text(3.6,5,'PI')
    text(6.6,-5,'PI')
    set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',[])
    t = title('Sensitivity Terms');
    set(t,'position',[4.4520 23.0000 -0.5000])
  end

  yyaxis right
  H = bar(8,co2_am_awm(i) .* 0.0423,'FaceColor','flat','BarWidth',.8,'EdgeColor','none');hold on
  H.CData = rgb('Orange');
  H = bar(7,co2_am_pi_awm(i) .* 0.0423,'FaceColor','flat','BarWidth',.8,'EdgeColor','none');hold on
  H.CData = rgb('OrangeRed');
  set(gca,'ytick',0:10:40)
  set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',[])
  ylim([0 40])
  ax = gca;
  ax.YColor = rgb('OrangeRed');
  ylabel('\muatm \circC^{-1}') 
  
  if i == 5
    set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',{'RF','AF','\delta\itp\rmCO_{2}'})
  else    
    set(gca,'xtick',[1.5 4.5 7.5],'xticklabel',[])
  end
end


%% Subplot Letters -----------------------------------

let = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)'};
for n = 1:length(let)
    subtightplot(5,3,n,[0.03 0.11],[0.05 .05],[0.11 0.08]);hold on
    
    yylim=get(gca,'YLim');xxlim=get(gca,'XLim');
    clear xr; xr = diff(xxlim);
    clear yr; yr = diff(yylim);

    if n == 3 | n == 6 | n == 9 | n == 12 | n == 15
      text(xxlim(1)-0.35.*xr,yylim(2)+.02.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left')
    elseif n == 1 | n == 4 | n == 7 | n == 10 | n == 13
      text(xxlim(1)-0.45.*xr,yylim(2)+.02.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left')
    else
      text(xxlim(1)-0.55.*xr,yylim(2)+.02.*yr,let(n),'VerticalAlignment','bottom','HorizontalAlignment','left')
    end
end


%% Save -----------------------------------

fig=gcf;
set(findall(fig,'-property','Fontname'),'fontname','Times','fontsize',9,'fontweight','normal');
exportgraphics(fig,[fig_path '/Figure_8.pdf'],'Resolution',300)
savefig([fig_path '/Figure_8'])


%% Clean up -----------------------------------

clear;close all;clc
