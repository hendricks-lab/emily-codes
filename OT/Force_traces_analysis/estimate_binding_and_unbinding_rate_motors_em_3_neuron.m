%adam g hendricks     27 June 2018   mcgill u.
%Emily modified November 8 2021
% estimate the unbinding rate from stationary optical trapping records
% following the method in Berger, Klumpp, and Lipowsky, Arxiv 2018
%20220402 corrected factor of 10 that the code was off in the ftm
%variable-Emily

clear all; close all; clc;

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(0,'defaulttextinterpreter','tex');
set(0,'DefaultFigureWindowStyle','docked');


addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

ab=[];

%Emily modifying for live cell data


for kch = 1:2
k_choose = kch;
if k_choose == 1
datdir = '/Volumes/Emily_2022/Omar_OT/U2OS/Good_force_traces/Force_mat/'; kc = 2;
col1=[0.5 0.5 0.5];
file_name = '*.mat'; 
elseif k_choose == 2
datdir = '/Volumes/Emily_2022/Omar_OT/U2OS/Good_force_traces/Force_mat/'; kc=2;
col1=[0.152941182255745 0.227450981736183 0.372549027204514];
file_name = '*.mat'; 
end

cmap = colormap(lines(5));

Fthreshold =0.5; %[pN] force threshold to define and event - 1.5, Emily modified from 0.7 to 0.5
Cstall_time=0.25; %[s] minimum stall plateau time Emily modified from 0.3 to 0.01
Cstall_vel=0; %[nm/s] minimum stall velocity (up)
Cstall_snapvel=20;% %[nm/s] minimum stall velocity (dwn) - 30 works really well Emily modified from 30 to 20
Cstall_time_bw_stall = 0.05; %[s] minimum time between stall events Emily tried 0.01 but 0.05 was better
zp = 0; %plot each trace? (0/1) Emily changed to 1 to see

dd = dir(fullfile(datdir,file_name)); %find all data files in datdir

force_bins = [-25:0.5:-Fthreshold,Fthreshold:0.5:25]; %Emily changed bins from -15 to 15 to -25 to 25
C_unbinding = zeros(size(force_bins)); %number of unbinding events in each force bin
C_points = zeros(size(force_bins)); %number of data points in each force bin
binding_rate = [];

Nf = numel(dd);
hw = waitbar(Nf,'Analyzing force traces ...');
for kf = 1:numel(dd);
    waitbar(kf/numel(dd),hw)
    clear datk
    datk = load(fullfile(datdir,dd(kf).name),'ftm','fpos','ffm','k_trap1');
%     datk=load([datdir,flist(kf).name]);
    
    %find number of 
    %find number of 
    datk.C_points = hist(datk.ffm,force_bins);
    
    %find unbinding events
    datk.ffm_smooth = sgolayfilt(datk.ffm, 1, 2001); %smoothing filter
    datk.fpos_smooth = sgolayfilt(datk.fpos, 1, 2001);
    
    [value1,index1]=max(abs(datk.ffm_smooth));
    
    datk.jstall=find(abs(datk.ffm_smooth)>Fthreshold);
    datk.ppm=datk.ffm_smooth./datk.k_trap1;
    if isempty(datk.jstall)==0;
    datk.jevent=find(diff(datk.ftm(datk.jstall)/10)>Cstall_time_bw_stall);
    datk.jevent=[0; datk.jevent; length(datk.jstall)];
    for kevent=2:length(datk.jevent);
        ji=datk.jstall(datk.jevent(kevent-1)+1);
        jf=datk.jstall(datk.jevent(kevent));
        
        
        %figure(2), hold on, plot(datk.ftm(ji:jf),Fthreshold.*ones(size(datk.ftm(ji:jf))),'c.')
        %hold on, plot(datk.ftm,datk.ffm,'Color','k');
        %hold on, plot(datk.ftm,datk.ffm_smooth,'Color','r','MarkerSize',20);
        
        stall_ti(kevent-1)=datk.ftm(ji)/10;
        stall_tf(kevent-1)=datk.ftm(jf)/10;
        [afmax, jmax] = max(abs(datk.ffm_smooth(ji:jf)));
        datk.stall_force_t(kevent-1) = datk.ftm(ji+jmax-1)/10;
        datk.stall_force(kevent-1)=datk.ffm_smooth(ji+jmax-1);%max(datk.ffm_smooth(ji:jf));
        
        datk.stall_time(kevent-1)=datk.ftm(jf)/10-datk.ftm(ji)/10; %plateau time
        
        %jvel=find(datk.ftm>(datk.ftm(ji)-Cstall_time) & datk.ftm<(datk.ftm(ji)+Cstall_time)); %%%%%%%%%%
        %pvel=polyfit(datk.ftm(jvel),datk.ppm(jvel),1);%%%%%%%%%%%%%%
        %datk.stall_vel(kevent-1)=pvel(1); % stall velocity %%%%%%%%%%%%
        jsnap=find(datk.ftm/10>(datk.ftm(jf)/10-0.05) & datk.ftm/10<(datk.ftm(jf)/10+0.05));%
        psnap=polyfit(datk.ftm(jsnap)/10,datk.fpos_smooth(jsnap),1);
        datk.stall_snapvel(kevent-1)=psnap(1); %snap back velocity
        if kevent > 2, %only calculate binding rate for events where there was a previous event in the trace
            datk.stall_binding_rate(kevent-1) = 1/(datk.ftm(ji)/10-datk.ftm(jf_previous_stall)/10);
        else
            datk.stall_binding_rate(kevent-1) = NaN;
        end
        jf_previous_stall = jf;
    end
    
    j_stall_keep = find(datk.stall_time > Cstall_time & abs(datk.stall_snapvel) > Cstall_snapvel);% & datk.stall_time < 20);% & abs(datk.stall_vel) > Cstall_vel);

    datk.C_unbinding = hist(datk.stall_force(j_stall_keep), force_bins);
    
    binding_rate = [binding_rate, [datk.stall_force(j_stall_keep); datk.stall_binding_rate(j_stall_keep)]];
    
    plus_binding_events=find(datk.stall_force(j_stall_keep)>0);
    minus_binding_events=find(datk.stall_force(j_stall_keep)<0);
    
    NumOfEvents_plus{kf} = numel(plus_binding_events);
    NumOfEvents_minus{kf} = numel(minus_binding_events);
    
    C_unbinding(kf,:) = datk.C_unbinding;
    C_points(kf,:) = datk.C_points;
    
    if zp==1,
       figure, plot(datk.ftm/10, datk.ffm), xlabel('time (s)'), ylabel('force (pN)')
       hold on, plot(datk.ftm/10, datk.ffm_smooth)
       hold on, plot(datk.stall_force_t(j_stall_keep)/10, datk.stall_force(j_stall_keep), 'ro', 'linewidth', 2)  
    end   
    end
end
close(hw)

dT = mean(diff(datk.ftm/10)); %[s] time b/w data points
Unbinding_Rate = sum(C_unbinding,1)./ (dT.*sum(C_points,1));

figure('Name',append(num2str(k_choose),'Unbinding force histogram plus'),'NumberTitle','off'), hold on, 
hau = bar(force_bins, sum(C_unbinding,1)); xlabel('Unbinding Force (pN)'), ylabel('Number of Events'), hau.FaceAlpha = 0.2; hau.FaceColor = cmap(kc,:);

figure('Name',append(num2str(k_choose),'Unbinding rate force dependence'),'NumberTitle','off'), hold on, 
hur = plot(force_bins, Unbinding_Rate, 'o','linewidth',2); xlabel('Force (pN)'), ylabel('Unbinding Rate (1/s)'); hur.Color = cmap(kc,:);

jb_plus = find(binding_rate(1,:) > 0 & isnan(binding_rate(2,:))==0);
jb_minus = find(binding_rate(1,:) < 0 & isnan(binding_rate(2,:)) == 0);

%% This is where you bin
% xb = 0.05:0.5:6; %original
xb = 0.05:0.3:20;

%mean_plus=mle(binding_rate(2,jb_plus)','distribution','exp');
%mean_minus=mle(binding_rate(2,jb_minus)','distribution','exp');
%prob_plus=pdf(mean_plus,xb);
%prob_minus=pdf(mean_minus,xb);

nb_plus = hist(binding_rate(2,jb_plus),xb);
nb_minus = hist(binding_rate(2,jb_minus),xb);

%fit_plus=fit(xb',nb_plus','exp1');
%fit_minus=fit(xb',nb_minus','exp1');

mean_plus=expfit(binding_rate(2,jb_plus));
mean_minus=expfit(binding_rate(2,jb_minus));

% figure, bar(xb,nb_plus./sum(nb_plus));
% [ft_plus,goodness_plus]=fit(xb',[nb_plus./sum(nb_plus)]','exp1');
% hold on, plot(ft_plus,xb,[nb_plus./sum(nb_plus)]);
% 
% figure, bar(xb,nb_minus./sum(nb_minus));
% [ft_minus,goodness_minus]=fit(xb',[nb_minus./sum(nb_minus)]','exp1');
% hold on, plot(ft_minus,xb,[nb_minus./sum(nb_minus)]);

Nb=sum(nb_plus) + sum(nb_minus);
figure('Name',append(num2str(k_choose),'Binding rate plus and minus histogram'),'NumberTitle','off'), hold on, 
subplot(1,2,1), hab = bar(xb,nb_plus./Nb); hab.FaceAlpha = 0.8; hab.FaceColor=col1;
hold on, subplot(1,2,1), ham = bar(xb, - nb_minus./Nb); ham.FaceAlpha = 0.8; ham.FaceColor = col1;
%hold on, subplot(1,2,1), plot(xb,prob_plus,'Color',col1,'LineWidth',2);
%hold on, subplot(1,2,1), plot(xb,prob_minus,'Color',col1,'LineWidth',2);
xlabel('Binding Rate (1/s)'), ylabel('Fraction of Events');
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
set(gca,'Ygrid','on');


ab(kch).avg_binding_rate_plus=mean_plus;
ab(kch).sem_binding_rate_plus=std(binding_rate(2,jb_plus))./numel(binding_rate(2,jb_plus));
ab(kch).avg_binding_rate_minus=mean_minus;
ab(kch).sem_binding_rate_minus=std(binding_rate(2,jb_minus))./numel(binding_rate(2,jb_minus));
ab(kch).binding_rate_plus=binding_rate(2,jb_plus);
ab(kch).binding_rate_minus=binding_rate(2,jb_minus);
ab(kch).unbinding_rate=Unbinding_Rate;
ab(kch).force_bins=force_bins;
ab(kch).NumOfEvents_plus=NumOfEvents_plus;
ab(kch).NumOfEvents_minus=NumOfEvents_minus;
%Emily added for plotting later
ab(kch).bins_binding=xb;
ab(kch).binned_binding_plus=nb_plus./Nb;
ab(kch).binned_binding_minus=-nb_minus./Nb;

clearvars -except ab
end

color_WT_plus=[0.5 0.5 0.5]; %Early
color_WT_minus=[0 0 0]; %Early
color_S421D_plus=[0 0.75 0.75]; %Early
color_S421D_minus=[0 0.4 0.4]; %Early

% 
% color_WT_plus=[0.2902 0.5569  0.8863]; %Late
% color_WT_minus=[0.1020    0.3490    0.6549]; %Late
% color_S421D_plus=[0.4314    0.6863    0.4745]; %Late
% color_S421D_minus=[0.2118    0.3843    0.2392]; %Late


%fname='Single-FL-kin1_';
%dir_save='G:\MAP7_v2\Single FL-Kinesin1\Processed data\Results\Binding rate\';
%save([dir_save,num2str(fname) '_Binding_rates.mat' ],'ab');

lw = 2;
kcol_black = color_WT_minus;
conf=0.95;  % 95% Confidence Interval

% Load all the values:
plus_ctrl=ab(1).binding_rate_plus;
plus_mean_ctrl=ab(1).avg_binding_rate_plus;
minus_ctrl=ab(1).binding_rate_minus;
minus_mean_ctrl=ab(1).avg_binding_rate_minus;
plus_S421D=ab(2).binding_rate_plus;
plus_mean_S421D=ab(2).avg_binding_rate_plus;
minus_S421D=ab(2).binding_rate_minus;
minus_mean_S421D=ab(2).avg_binding_rate_minus;

%Emily adding unbinding rate plots (Apr 3/2022)
% unbinding_plus_WT=ab(1).unbinding_rate(find(ab(1).force_bins>0));
% unbinding_minus_WT=ab(1).unbinding_rate(find(ab(1).force_bins<0));
% unbinding_plus_S421D=ab(2).unbinding_rate(find(ab(2).force_bins>0));
% unbinding_minus_S421D=ab(2).unbinding_rate(find(ab(2).force_bins<0));
% [plus_cdf_WT, x_cdf_plus_WT]=ecdf(unbinding_plus_WT);
% [minus_cdf_WT, x_cdf_minus_WT]=ecdf(unbinding_minus_WT);
% [plus_cdf_S421D, x_cdf_plus_S421D]=ecdf(unbinding_plus_S421D);
% [minus_cdf_S421D, x_cdf_minus_S421D]=ecdf(unbinding_minus_S421D);
% 
% 
% jj_plus_WT = 0.1.*randn(size(unbinding_plus_WT));
% jj_minus_WT = 2+0.1.*randn(size(unbinding_minus_WT));
% jj_plus_S421D = 1+0.1.*randn(size(unbinding_plus_S421D));
% jj_minus_S421D = 3+0.1.*randn(size(unbinding_minus_S421D));
% 
% %assign these variables and add errorbars to scatterplot
% sem_plus_WT=nanstd(unbinding_plus_WT)./sqrt(numel(unbinding_plus_WT));
% sem_minus_WT=nanstd(unbinding_minus_WT)./sqrt(numel(unbinding_minus_WT));
% mean_plus_WT=nanmean(unbinding_plus_WT);
% mean_minus_WT=nanmean(unbinding_minus_WT);
% sem_plus_S421D=nanstd(unbinding_plus_S421D)./sqrt(numel(unbinding_plus_S421D));
% sem_minus_S421D=nanstd(unbinding_minus_S421D)./sqrt(numel(unbinding_minus_S421D));
% mean_plus_S421D=nanmean(unbinding_plus_S421D);
% mean_minus_S421D=nanmean(unbinding_minus_S421D);
% 
% 
% figure('Name','Unbinding rate scatter both','NumberTitle','off'), hold on, 
% plot(jj_plus_WT,unbinding_plus_WT, 'o', 'Color', color_WT_plus); %WT
% hold on, plot(jj_minus_WT,unbinding_minus_WT, 'o', 'Color', color_WT_minus);%WT
% hold on, errorbar(0,mean_plus_WT,sem_plus_WT,'CapSize',8,'LineWidth',2,'Color','k')
% hold on, errorbar(2,mean_minus_WT,sem_minus_WT,'CapSize',8,'LineWidth',2,'Color','k')
% plot(jj_plus_S421D,unbinding_plus_S421D, 'o', 'Color', color_S421D_plus);%S421D
% hold on, plot(jj_minus_S421D,unbinding_minus_S421D, 'o', 'Color', color_S421D_minus);%S421D
% errorbar(1,mean_plus_S421D,sem_plus_S421D,'CapSize',8,'LineWidth',2,'Color','k')
% errorbar(3,mean_minus_S421D,sem_minus_S421D,'CapSize',8,'LineWidth',2,'Color','k')
% set(gca, 'yscale', 'log');
% set(gca,'xticklabel',{'WT Plus','S421D Plus','WT Minus', 'S421D Minus'});
% set(gca,'xlim',[-0.5 3.5]);
% set(gca, 'xtick', [0 1 2 3]);
% ylabel('Unbinding Rate (1/s)');
%   set(gca,'FontSize',20);
%   set(gca, 'FontWeight','bold');
%   set(gca, 'LineWidth',2);
%   pbaspect([1 1 1]);
%   box on
% 
% figure('Name','Unbinding rate cdf both','NumberTitle','off'), hold on, 
% plot(plus_cdf_WT, x_cdf_plus_WT,'-','linewidth',2,'Color',color_WT_plus);
% hold on, plot(minus_cdf_WT,x_cdf_minus_WT,'-','linewidth',2,'Color',color_WT_minus); 
% plot(plus_cdf_S421D, x_cdf_plus_S421D,'-','linewidth',2,'Color',color_S421D_plus); %S421D
% hold on, plot(minus_cdf_S421D,x_cdf_minus_S421D,'-','linewidth',2,'Color',color_S421D_minus); %S421D
% set(gca, 'yscale', 'log');
%   ylabel('Unbinding Rate (1/s)');
%   set(gca,'FontSize',20);
%   set(gca, 'FontWeight','bold');
%   set(gca, 'LineWidth',2);
%   pbaspect([1 1 1]);
%   box on
%End of Emily adding unbinding rate plots
% Calculate the Confidence Intervals:
% ci_plus_ctrl=ConfidenceInterval(plus_ctrl,conf);
% ci_plus_S421D=ConfidenceInterval(plus_S421D,conf);
% ci_minus_ctrl=ConfidenceInterval(minus_ctrl,conf);
% ci_minus_S421D=ConfidenceInterval(minus_S421D,conf);

figure('Name','Box plot binding rate plus and minus both k_chooses','NumberTitle','off'), hold on, 
subplot(1,2,1), hold on, 
h1=notBoxPlot(plus_ctrl,1-0.8,'style','line');
set(h1.data,'MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
h2=notBoxPlot(plus_S421D,1+0.8,'style','line');
set(h2.data,'MarkerFaceColor',[0.152941182255745 0.227450981736183 0.372549027204514],'MarkerEdgeColor',[0.152941182255745 0.227450981736183 0.372549027204514]);
set(gca,'Ygrid','on');
set(gca,'LineWidth',2);
set(gca,'FontSize',20);
xlabel(''); ylabel('Plus Direction Binding rate (s^{-1})');
xlim([-1 3]);
set(gca,'xticklabel',{[]}) 
%xticks(0.7,1,3);
%xticklabels({'WT','S421D'});
box on;

subplot(1,2,2), hold on, 
h3=notBoxPlot(minus_ctrl,1-0.8,'style','line');
set(h3.data,'MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
h4=notBoxPlot(minus_S421D,1+0.8,'style','line');
set(h4.data,'MarkerFaceColor',[0.152941182255745 0.227450981736183 0.372549027204514],'MarkerEdgeColor',[0.152941182255745 0.227450981736183 0.372549027204514]);
set(gca,'Ygrid','on');
set(gca,'LineWidth',2);
set(gca,'FontSize',20);
xlabel(''); ylabel('Minus Direction Binding Rate (s^{-1})');
xlim([-1 3]);
set(gca,'xticklabel',{[]}) 
%xticks(0.7,1,3);
%xticklabels({'WT','S421D'});
box on;

figure('Name','Binding box plot plus and minus both k_chooses','NumberTitle','off'), hold on,  
h1=notBoxPlot(plus_ctrl,1-0.8,'style','line');
set(h1.data,'MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
h2=notBoxPlot(plus_S421D,0.8,'style','line');
set(h2.data,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor',color_S421D_plus);
h3=notBoxPlot(-minus_ctrl,1-0.8,'style','line');
set(h3.data,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor',color_WT_minus);
hold on,
h4=notBoxPlot(-minus_S421D,0.8,'style','line');
set(h4.data,'MarkerFaceColor',color_S421D_minus,'MarkerEdgeColor',color_S421D_minus);
hold on,
plot([-0.2;1.2],[0;0],'k-');
set(gca,'Ygrid','off');
set(gca,'FontSize',20);
xlabel({'WT                          S421D','Late Endosomes'}); 
ylabel('Binding rate (s^{-1})');
xlim([-0.2 1.2]);
set(gca,'xticklabel',{[]});
publication_fig(0,0,1);

figure('Name','Binding box plot plus and minus both k_chooses','NumberTitle','off'), hold on,  
h1=notBoxPlot(plus_ctrl,1-0.8,'style','line');
set(h1.data,'MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
h2=notBoxPlot(plus_S421D,0.8,'style','line');
set(h2.data,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor',color_S421D_plus);
h3=notBoxPlot(-minus_ctrl,1-0.8,'style','line');
set(h3.data,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor',color_WT_minus);
hold on,
h4=notBoxPlot(-minus_S421D,0.8,'style','line');
set(h4.data,'MarkerFaceColor',color_S421D_minus,'MarkerEdgeColor',color_S421D_minus);
hold on,
plot([-0.2;1.2],[0;0],'k-');
set(gca,'Ygrid','off');
set(gca,'FontSize',20);
xlabel({'WT                          S421D','Early Endosomes'}); 
ylabel('Binding rate (s^{-1})');
xlim([-0.2 1.2]);
set(gca,'xticklabel',{[]});
publication_fig(0,0,1);

figure('Name','Binding rate plus and minus histogram both conditions','NumberTitle','off'), hold on, 
hab = bar(ab(1).bins_binding,ab(1).binned_binding_plus); hab.FaceAlpha = 0.5; hab.FaceColor=color_WT_minus; hab.EdgeColor=color_WT_minus;
hold on, ham = bar(ab(1).bins_binding, ab(1).binned_binding_minus); ham.FaceAlpha = 0.5; ham.FaceColor = color_WT_minus; ham.EdgeColor=color_WT_minus;
habc2 = bar(ab(2).bins_binding,ab(2).binned_binding_plus); habc2.FaceAlpha = 0.5; habc2.FaceColor=color_S421D_minus; habc2.EdgeColor=color_S421D_minus;
hold on, hamc2 = bar(ab(2).bins_binding, ab(2).binned_binding_minus); hamc2.FaceAlpha = 0.5; hamc2.FaceColor = color_S421D_minus; hamc2.EdgeColor=color_S421D_minus;
xlabel({'Binding Rate (1/s)','Early Endosomes'}), ylabel('Fraction of Events');
xlim([0 20]);
publication_fig(0,0,1);


figure('Name','Scatter plus and minus binding rate both k_chooses','NumberTitle','off'), hold on, 
subplot(1,2,1), hold on, 
sem_pl_ctrl=std(plus_ctrl)./sqrt(numel(plus_ctrl));
sem_pl_S421D=std(plus_S421D)./sqrt(numel(plus_S421D));
plot(1-0.8,plus_mean_ctrl,'o','MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
errorbar(1-0.8,plus_mean_ctrl,sem_pl_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
hold on,
plot(1+0.8,plus_mean_S421D,'o','MarkerEdgeColor',color_S421D_minus);
hold on,
errorbar(1+0.8,plus_mean_S421D,sem_pl_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);

subplot(1,2,2), hold on, 
sem_mi_ctrl=std(minus_ctrl)./sqrt(numel(minus_ctrl));
sem_mi_S421D=std(minus_S421D)./sqrt(numel(minus_S421D));
plot(1-0.8,minus_mean_ctrl,'o','MarkerFaceColor',color_WT_plus,'MarkerEdgeColor',color_WT_plus);
hold on,
errorbar(1-0.8,minus_mean_ctrl,sem_mi_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
hold on,
plot(1+0.8,minus_mean_S421D,'o','MarkerEdgeColor',color_S421D_minus);
hold on,
errorbar(1+0.8,minus_mean_S421D,sem_mi_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);

figure('Name','Line and SEM binding rate plus and minus both condiitons','NumberTitle','off'), hold on, 
line([ab(1).avg_binding_rate_plus,ab(1).avg_binding_rate_plus],[0,1]);
hold on, 
line([ab(2).avg_binding_rate_plus,ab(2).avg_binding_rate_plus],[0,1]);

[kstest_plus_WT_vs_S421D p_kstest_plus_WT_vs_S421D]=kstest2(plus_ctrl,plus_S421D)
[kstest_minus_WT_vs_S421D p_kstest_minus_WT_vs_S421D]=kstest2(minus_ctrl,minus_S421D)
[kstest_plus_vs_minus_WT p_kstest_plus_vs_minus_WT]=kstest2(plus_ctrl,minus_ctrl)
[kstest_plus_vs_minus_S421D p_kstest_plus_vs_minus_S421D]=kstest2(plus_S421D,minus_S421D)

datplus = {plus_ctrl, plus_S421D};
datminus = {minus_ctrl, minus_S421D};


%Calculating the cumulative distribution function for each of the binding
%rates in both plus and minus directions

jj=cell(2,1);
jj{1} = 0.1.*randn(size(datplus{1}));
jj{2} = 1+0.1.*randn(size(datplus{2}));

rr = logspace(-2,1.5,10);
nr1 = hist(datplus{1},rr);
nr2 = hist(datplus{2},rr);

[c1, x1] = ecdf(datplus{1});
[c2, x2] = ecdf(datplus{2});

sd_plus_ctrl=std(plus_ctrl);
sd_minus_ctrl=std(minus_ctrl);
sd_plus_S421D=std(plus_S421D);
sd_minus_S421D=std(minus_S421D);

figure('Name','Plus binding rate with CDF','NumberTitle','off'), hold on, 
  subplot(121), hold on
  h11=plot(jj{1},datplus{1},'o','MarkerEdgeColor',color_WT_minus);
%   errorbar(0,plus_mean_ctrl,sd_plus_ctrl,'CapSize',8,'LineWidth',lw,'Color',color_WT_minus);
  errorbar(0,plus_mean_ctrl,sem_pl_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  h12=plot(jj{2},datplus{2},'o','MarkerEdgeColor',color_S421D_minus);
%   errorbar(1,plus_mean_S421D,sd_plus_S421D,'CapSize',8,'LineWidth',lw,'Color',color_WT_minus);
  errorbar(1,plus_mean_S421D,sem_pl_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  a1=gca;
  set(a1,'YScale','log');
  set(a1,'ylim',[3e-2, 3.1e1]);
  set(a1,'xtick',[0 1]);
  set(a1,'xticklabel',{'CTRL','S421D'});
  ylabel('Binding Rate (1/s)');
  set(a1,'FontSize',20);
  set(a1, 'FontWeight','bold');
  set(a1, 'LineWidth',2);
  pbaspect([1 1 1]);
  box on

  subplot(122), hold on
  %h21=plot(cumsum(nr1./sum(nr1)),rr,'-','linewidth',2)
  %h22=plot(cumnr2./sum(nr2),rr,'-','linewidth',2)
  h21 = plot(c1, x1, '-','linewidth',2,'Color',color_WT_minus);
  h22 = plot(c2, x2, '-','linewidth',2,'Color',color_S421D_minus);
  a2=gca;
  set(a2,'yscale','log')
  set(a2,'ylim',[3e-2, 3.1e1])
  set(a2,'FontSize',20);
  set(a2, 'FontWeight','bold');
  set(a2, 'LineWidth',2);
  pbaspect([1 1 1]);
  xlabel('CDF')
  box on
  
jjminus=cell(2,1);
jjminus{1} = 2+0.1.*randn(size(datminus{1}));
jjminus{2} = 3+0.1.*randn(size(datminus{2}));

rrminus = logspace(-2,1.5,10);
nr1minus = hist(datminus{1},rrminus);
nr2minus = hist(datminus{2},rrminus);

[c1minus, x1minus] = ecdf(datminus{1});
[c2minus, x2minus] = ecdf(datminus{2});
  
  figure('Name','Binding rate with CDF','NumberTitle','off'), hold on, %Emily putting all binding rates together
  subplot(121), hold on
  plot(jj{1},datplus{1},'o','MarkerEdgeColor',color_WT_plus);
  errorbar(0,plus_mean_ctrl,sem_pl_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  plot(jjminus{1},datminus{1},'o','MarkerEdgeColor',color_WT_minus);
  errorbar(2,minus_mean_ctrl,sem_mi_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  plot(jj{2},datplus{2},'o','MarkerEdgeColor',color_S421D_plus);
  plot(jjminus{2},datminus{2},'o','MarkerEdgeColor',color_S421D_minus);
  errorbar(1,plus_mean_S421D,sem_pl_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  errorbar(3,minus_mean_S421D,sem_mi_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  a1=gca;
  set(a1,'YScale','log');
  set(a1,'ylim',[3e-2, 3.1e1]);
  set(gca,'xticklabel',{'WT+','WT-','S421D+', 'S421D-'});
  set(gca,'xlim',[-0.5 3.5]);
  set(gca, 'xtick', [0 1 2 3]);
  ylabel('Binding Rate (1/s)');
  set(a1,'FontSize',20);
%   set(a1, 'FontWeight','bold');
  set(a1, 'LineWidth',2);
%   pbaspect([1 1 1]);
  box on

  subplot(122), hold on
  %h21=plot(cumsum(nr1./sum(nr1)),rr,'-','linewidth',2)
  %h22=plot(cumnr2./sum(nr2),rr,'-','linewidth',2)
  plot(c1, x1, '-','linewidth',2,'Color',color_WT_plus);
  plot(c2, x2, '-','linewidth',2,'Color',color_S421D_plus);
  plot(c1minus, x1minus, '-','linewidth',2,'Color',color_WT_minus);
  plot(c2minus, x2minus, '-','linewidth',2,'Color',color_S421D_minus);
  a2=gca;
  set(a2,'yscale','log')
  set(a2,'ylim',[3e-2, 3.1e1])
  set(a2,'FontSize',20);
%   set(a2, 'FontWeight','bold');
  set(a2, 'LineWidth',2);
%   pbaspect([1 1 1]);
  xlabel('CDF')
  box on
  
  %Early plots
  figure('Name','Binding rate scatter','NumberTitle','off'), hold on, %Emily putting all binding rates together
  plot(jj{1},datplus{1},'o','MarkerEdgeColor',color_WT_plus);
  errorbar(0,plus_mean_ctrl,sem_pl_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  plot(jjminus{1},datminus{1},'o','MarkerEdgeColor',color_WT_minus);
  errorbar(2,minus_mean_ctrl,sem_mi_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  plot(jj{2},datplus{2},'o','MarkerEdgeColor',color_S421D_plus);
  plot(jjminus{2},datminus{2},'o','MarkerEdgeColor',color_S421D_minus);
  errorbar(1,plus_mean_S421D,sem_pl_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  errorbar(3,minus_mean_S421D,sem_mi_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  publication_fig(0,0,1)
  a1=gca;
%   set(a1,'YScale','log');
%   set(a1,'ylim',[3e-2, 3.1e1]);
  set(gca,'xticklabel',{'WT Plus','WT Minus','S421D Plus', 'S421D Minus'});
  set(gca,'xlim',[-0.5 3.5]);
  set(gca, 'xtick', [0 1 2 3]);
  ylabel('Binding Rate (1/s)');
  set(a1,'FontSize',20);
  set(a1, 'LineWidth',2);
  pbaspect([1 1 1]);
  box on

  figure('Name','Binding rate cdf','NumberTitle','off'), hold on, %Emily plotting cdf separate to have square plots
  plot(c1, x1, '-','linewidth',2,'Color',color_WT_plus);
  plot(c2, x2, '-','linewidth',2,'Color',color_S421D_plus);
  plot(c1minus, x1minus, '-','linewidth',2,'Color',color_WT_minus);
  plot(c2minus, x2minus, '-','linewidth',2,'Color',color_S421D_minus);
  a2=gca;
%   set(a2,'yscale','log')
%   set(a2,'ylim',[3e-2, 3.1e1])
  set(a2,'FontSize',20);
  ylabel('Binding Rate (1/s)');
  set(a2, 'LineWidth',2);
  pbaspect([1 1 1]);
  xlabel('CDF')
  box on
  
  %Late plots
%   figure('Name','Binding rate scatter','NumberTitle','off'), hold on, %Emily putting all binding rates together
%   plot(jj{1},datplus{1},'o','MarkerEdgeColor',[0.2902 0.5569  0.8863]);
%   errorbar(0,plus_mean_ctrl,sem_pl_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
%   plot(jjminus{1},datminus{1},'o','MarkerEdgeColor',[0.1020    0.3490    0.6549]);
%   errorbar(2,minus_mean_ctrl,sem_mi_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
%   plot(jj{2},datplus{2},'o','MarkerEdgeColor',[0.4314    0.6863    0.4745]);
%   plot(jjminus{2},datminus{2},'o','MarkerEdgeColor',[0.2118    0.3843    0.2392]);
%   errorbar(1,plus_mean_S421D,sem_pl_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
%   errorbar(3,minus_mean_S421D,sem_mi_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
%   a1=gca;
% %   set(a1,'YScale','log');
% %   set(a1,'ylim',[3e-2, 3.1e1]);
%   set(gca,'xticklabel',{'WT Plus','WT Minus','S421D Plus', 'S421D Minus'});
%   set(gca,'xlim',[-0.5 3.5]);
%   set(gca, 'xtick', [0 1 2 3]);
%   ylabel('Binding Rate (1/s)');
%   set(a1,'FontSize',20);
%   set(a1, 'LineWidth',2);
%   pbaspect([1 1 1]);
%   box on
% 
%   figure('Name','Binding rate cdf','NumberTitle','off'), hold on, %Emily plotting cdf separate to have square plots
%   plot(c1, x1, '-','linewidth',2,'Color',[0.2902 0.5569  0.8863]);
%   plot(c2, x2, '-','linewidth',2,'Color',[0.4314    0.6863    0.4745]);
%   plot(c1minus, x1minus, '-','linewidth',2,'Color',[ 0.1020    0.3490    0.6549]);
%   plot(c2minus, x2minus, '-','linewidth',2,'Color',[0.2118    0.3843    0.2392]);
%   a2=gca;
% %   set(a2,'yscale','log')
% %   set(a2,'ylim',[3e-2, 3.1e1])
%   set(a2,'FontSize',20);
%   ylabel('Binding Rate (1/s)');
%   set(a2, 'LineWidth',2);
%   pbaspect([1 1 1]);
%   xlabel('CDF')
%   box on
  
  
% jjminus=cell(2,1);
% jjminus{1} = 0.1.*randn(size(datminus{1}));
% jjminus{2} = 1+0.1.*randn(size(datminus{2}));
% 
% rrminus = logspace(-2,1.5,10);
% nr1minus = hist(datminus{1},rrminus);
% nr2minus = hist(datminus{2},rrminus);
% 
% [c1minus, x1minus] = ecdf(datminus{1});
% [c2minus, x2minus] = ecdf(datminus{2});

%Emily look up how to make the error bar 2 sided
figure('Name','Minus binding rate with CDF','NumberTitle','off'), hold on,  
  subplot(121), hold on
  h11=plot(jjminus{1},datminus{1},'o','MarkerEdgeColor',color_WT_minus);
%   h111=errorbar(0,minus_mean_ctrl,sd_minus_ctrl,'CapSize',8,'LineWidth',lw,'Color',color_WT_minus);
  h112=errorbar(0,minus_mean_ctrl,sem_mi_ctrl,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  h12=plot(jjminus{2},datminus{2},'o','MarkerEdgeColor',color_S421D_minus);
%   h121=errorbar(1,minus_mean_S421D,sd_minus_S421D,'CapSize',8,'LineWidth',lw,'Color',color_WT_minus);
  h122=errorbar(1,minus_mean_S421D,sem_mi_S421D,'CapSize',8,'LineWidth',lw,'Color',kcol_black);
  a1=gca;
  set(a1,'YScale','log');
  set(a1,'ylim',[3e-2, 3.1e1]);
  set(a1,'xtick',[0 1]);
  set(a1,'xticklabel',{'CTRL','S421D'});
  ylabel('Binding Rate (1/s)');
  set(a1,'FontSize',20);
  set(a1, 'FontWeight','bold');
  set(a1, 'LineWidth',2);
  pbaspect([1 1 1]);
  box on

  subplot(122), hold on
  %h21=plot(cumsum(nr1./sum(nr1)),rr,'-','linewidth',2)
  %h22=plot(cumnr2./sum(nr2),rr,'-','linewidth',2)
  h21 = plot(c1minus, x1minus, '-','linewidth',2,'Color',color_WT_minus);
  h22 = plot(c2minus, x2minus, '-','linewidth',2,'Color',color_S421D_minus);
  a2=gca;
  set(a2,'yscale','log');
  set(a2,'ylim',[3e-2, 3.1e1]);
  xlabel('CDF');
  set(a2,'FontSize',20);
  set(a2, 'FontWeight','bold');
  set(a2, 'LineWidth',2);
  pbaspect([1 1 1]);
  box on
   
% 
fileprefix='20240419_U2OS_stall_0pt1_binding_';
% 
tempdir = '/Volumes/Emily_2022/Omar_OT/U2OS/Good_force_traces/';   % Your destination folder
FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end
