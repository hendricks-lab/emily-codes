%% Analyze velocties and run length of processive events:
% Code #3 (final) in the per trajectory analysis.
%Steps to run this code:
% 1.Update the addpaths (lines 11-16) and the directories (lines 33-110), and your DT.
% 2.Update the directory for saving figures (if you want to save them),
% lines (618-620).
% 3.Run the code!

clear all; close all; clc; 
cd(strtok(userpath, pathsep));
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

set(0,'defaulttextinterpreter','tex');
set(0,'DefaultFigureWindowStyle','docked');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

%Setting the colour code for each condition (because it's so annoying to do
%in each plot)
% polyQ length
cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_45Q=cmap(20,:);%[0.75 0 0];
colour_65Q=cmap(25,:);%[0.5 0 0];
colour_81Q=cmap(30,:);%[0.25 0 0];\
htt_30Q="30Q";
htt_45Q="45Q";
htt_65Q="65Q";
htt_81Q="81Q";

%stress condition
% cmap = colormap(hot(100));
% colour_30Q=[0.5 0.5 0.5];
% colour_65Q=cmap(30,:);%[0.5 0 0];
% cmap2= colormap(spring(100));
% colour_45Q= [1.0000 0.4444 0.5556];
% colour_81Q=cmap2(30,:);%[0.25 0 0];
% htt_30Q="30Q";
% htt_45Q="30Q+IFg";
% htt_65Q="81Q";
% htt_81Q="81Q+IFg";

for k_choose = 1:4

if k_choose == 1 %30Q 
%     cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_pos/');
    cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_pos/');
    DT=0.12;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    fraction_plus_proc_time=[]; %Emily added
elseif k_choose == 2 %45Q 
%     cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_pos/');
%     cd('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_pos/');
    cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_pos/')
    DT=0.12;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    fraction_plus_proc_time=[]; %Emily added
elseif k_choose == 3 %65Q
%     cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_pos/');
%     cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_pos/');
    cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_pos/');
    DT=0.12;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    fraction_plus_proc_time=[]; %Emily added
elseif k_choose == 4 %81Q
%         cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_pos/');
%     cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_pos/');
%     cd('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_pos/');
    cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_pos/');
    DT=0.12;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    fraction_plus_proc_time=[]; %Emily added
end
min_lr=3;
k_choose
load('position.mat');

% end
for i=1:numel(position)
    pos1=position{i}; 
    time=[1:numel(pos1)]'.*DT;
    diff_coeff=Diffusion_coefficient(i);
    event_dur=time(end)-time(1);
    pos=pos1;%smooth(pos1,10,'sgolay',1);
    res=analyze_run_length_reversals_v5_per_cell_bias(time,pos,min_lr,0,k_choose);
    
    rev_rate=[rev_rate,(numel(res.proc_run)./(event_dur))];%[rev_rate,res.reversal_rate];
    proc_run=[proc_run;res.proc_run];
    proc_time=[proc_time;res.proc_time];
    diff_run=[diff_run;res.diff_run];
    diff_time=[diff_time;res.diff_time];
    event_dur_v2=[event_dur_v2,event_dur];
    
    fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time(find(res.proc_run>0)))/sum(res.proc_time))]; %Emily added to find plus end runs
    
    avg_lr_per_traj1p=[avg_lr_per_traj1p;res.avg_lrplus];
    avg_lr_per_traj1m=[avg_lr_per_traj1m;res.avg_lrminus];
    avg_vel_per_traj1p=[avg_vel_per_traj1p;res.avg_velplus];
    avg_vel_per_traj1m=[avg_vel_per_traj1m;res.avg_velminus];
    
    %pe1=func_peclet_number(res.run_bw_rev,res.time_bw_rev,Diffusion_coefficient(i));
    %pe=[pe,pe1];
%     dc=[dc,diff_coeff];
    all_pos=[all_pos;pos];
    
    %[deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d (xk1{i}, yk1{i}, 0.2, k_choose);
    %delt=[delt,deltat];
    %msd_pts=[msd_pts,msdpts];

    clear res pe1 deltat msdpts sem log_deltat log_msdpts alpha_1 Diff_Coef
end
disp(i)

proc_runs{k_choose}=proc_run;
proc_times{k_choose}=proc_time;
rev_rates{k_choose}=rev_rate;
diff_runs{k_choose}=diff_run;
diff_times{k_choose}=diff_time;
avg_lr_per_trajp{k_choose}=avg_lr_per_traj1p;
avg_vel_per_trajp{k_choose}=avg_vel_per_traj1p;
avg_lr_per_trajm{k_choose}=avg_lr_per_traj1m;
avg_vel_per_trajm{k_choose}=avg_vel_per_traj1m;
%peclet{k_choose}=pe;
diff_coeff_fn{k_choose}=dc;
pos_all{k_choose}=all_pos;
all_delt{k_choose}=delt;
all_msd{k_choose}=msd_pts;
% xk_all{k_choose}=xk1;
% yk_all{k_choose}=yk1;
posit{k_choose}=position;
duration_traj{k_choose}=event_dur_v2;

plus_proc_times{k_choose}=fraction_plus_proc_time; %Emily Added

clear position xk_traj yk_traj
clear pos time 
proc_run=[];
proc_time=[];
rev_rate=[];
diff_run=[];
diff_times=[];
pe=[];
dc=[];
all_pos=[];
event_dur_v2=[];
avg_lr_per_traj1p=[];
avg_lr_per_traj1m=[];
avg_vel_per_traj1p=[];
avg_vel_per_traj1m=[];
end

%% Analyze Reversal Rate:
% xk_rev=-1:0.01:1;
% rev_rate_18Q=rev_rates{1};
% rev_rate_45Q=rev_rates{2};
% hist_xkrev_18Q=hist(rev_rate_18Q,xk_rev);
% hist_xkrev_45Q=hist(rev_rate_45Q,xk_rev);
% figure(1), hold on, 
% bar(xk_rev,hist_xkrev_45Q./sum(hist_xkrev_45Q),'BarWidth',1,'FaceColor','g','EdgeColor','g');
% hold on, 
% stairs(xk_rev,hist_xkrev_18Q./sum(hist_xkrev_18Q),'Color','k','LineWidth',2);
% xlabel('Reversal Rate (# of reversals/s)'); ylabel('Fraction of Events'); 
% publication_fig(0,0,1);
% % ./sum(hist_xkrev_18Q+hist_xkrev_45Q)

%% Analyze Processive Runs:
% proc_runs_18Q=proc_runs{1};
% proc_times_18Q=proc_times{1};
proc_runs_30Q=proc_runs{1};
proc_times_30Q=proc_times{1};
proc_runs_45Q=proc_runs{2};
proc_times_45Q=proc_times{2};
proc_runs_65Q=proc_runs{3};
proc_times_65Q=proc_times{3};
proc_runs_81Q=proc_runs{4};
proc_times_81Q=proc_times{4};


% proc_plus_run_18Q=proc_runs_18Q(proc_runs_18Q>0);
% proc_plus_time_18Q=proc_times_18Q(proc_runs_18Q>0);
% proc_minus_run_18Q=proc_runs_18Q(proc_runs_18Q<0);
% proc_minus_time_18Q=proc_times_18Q(proc_runs_18Q<0);

proc_plus_run_30Q=proc_runs_30Q(proc_runs_30Q>0);
proc_plus_time_30Q=proc_times_30Q(proc_runs_30Q>0);
proc_minus_run_30Q=proc_runs_30Q(proc_runs_30Q<0);
proc_minus_time_30Q=proc_times_30Q(proc_runs_30Q<0);

proc_plus_run_45Q=proc_runs_45Q(proc_runs_45Q>0);
proc_plus_time_45Q=proc_times_45Q(proc_runs_45Q>0);
proc_minus_run_45Q=proc_runs_45Q(proc_runs_45Q<0);
proc_minus_time_45Q=proc_times_45Q(proc_runs_45Q<0);

proc_plus_run_65Q=proc_runs_65Q(proc_runs_65Q>0);
proc_plus_time_65Q=proc_times_65Q(proc_runs_65Q>0);
proc_minus_run_65Q=proc_runs_65Q(proc_runs_65Q<0);
proc_minus_time_65Q=proc_times_65Q(proc_runs_65Q<0);

proc_plus_run_81Q=proc_runs_81Q(proc_runs_81Q>0);
proc_plus_time_81Q=proc_times_81Q(proc_runs_81Q>0);
proc_minus_run_81Q=proc_runs_81Q(proc_runs_81Q<0);
proc_minus_time_81Q=proc_times_81Q(proc_runs_81Q<0);

% total_time_18Q=sum(proc_plus_time_18Q)+sum(proc_minus_time_18Q);
total_time_30Q=sum(proc_plus_time_30Q)+sum(proc_minus_time_30Q);
total_time_45Q=sum(proc_plus_time_45Q)+sum(proc_minus_time_45Q);
total_time_65Q=sum(proc_plus_time_65Q)+sum(proc_minus_time_65Q);
total_time_81Q=sum(proc_plus_time_81Q)+sum(proc_minus_time_81Q);


% fraction_plus_18Q=sum(proc_plus_time_18Q)./total_time_18Q;
% fraction_minus_18Q=sum(proc_minus_time_18Q)./total_time_18Q;
fraction_plus_30Q=sum(proc_plus_time_30Q)./total_time_30Q;
fraction_minus_30Q=sum(proc_minus_time_30Q)./total_time_30Q;
fraction_plus_45Q=sum(proc_plus_time_45Q)./total_time_45Q;
fraction_minus_45Q=sum(proc_minus_time_45Q)./total_time_45Q;
fraction_plus_65Q=sum(proc_plus_time_65Q)./total_time_65Q;
fraction_minus_65Q=sum(proc_minus_time_65Q)./total_time_65Q;
fraction_plus_81Q=sum(proc_plus_time_81Q)./total_time_81Q;
fraction_minus_81Q=sum(proc_minus_time_81Q)./total_time_81Q;


figure('Name','Directional Bias','NumberTitle','off'), hold on,
% bar(1.2,fraction_plus_18Q,'BarWidth',0.2,'FaceColor',colour_18Q,'EdgeColor',colour_18Q,'Facealpha',1);
% hold on, 
% bar(1.2,-1.*fraction_minus_18Q,'BarWidth',0.2,'FaceColor',colour_18Q,'EdgeColor',colour_18Q,'Facealpha',0.3);
% hold on, 
bar(1.4,fraction_plus_30Q,'BarWidth',0.2,'FaceColor',colour_30Q,'EdgeColor',colour_30Q,'Facealpha',1);
hold on, 
bar(1.4,-1.*fraction_minus_30Q,'BarWidth',0.2,'FaceColor', colour_30Q,'EdgeColor',colour_30Q,'Facealpha',0.3);
hold on, 
bar(1.6,fraction_plus_45Q,'BarWidth',0.2,'FaceColor',colour_45Q,'EdgeColor',colour_45Q,'Facealpha',1);
hold on, 
bar(1.6,-1.*fraction_minus_45Q,'BarWidth',0.2,'FaceColor',colour_45Q,'EdgeColor',colour_45Q,'Facealpha',0.3);
% hold on, 
bar(1.8,fraction_plus_65Q,'BarWidth',0.2,'FaceColor',colour_65Q,'EdgeColor',colour_65Q,'Facealpha',1);
hold on, 
bar(1.8,-1.*fraction_minus_65Q,'BarWidth',0.2,'FaceColor',colour_65Q,'EdgeColor',colour_65Q,'Facealpha',0.3);
hold on, 
bar(2,fraction_plus_81Q,'BarWidth',0.2,'FaceColor',colour_81Q,'EdgeColor',colour_81Q,'Facealpha',1);
hold on, 
bar(2,-1.*fraction_minus_81Q,'BarWidth',0.2,'FaceColor',colour_81Q,'EdgeColor',colour_81Q,'Facealpha',0.3);
ylabel('Fraction of time of Processive runs'); 
xlim([1 2.25]);
set(gca,'ygrid','on');
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2,'XColor',[0 0 0],...
    'XTick',[1.25 1.5 1.75 2 2.25],'XTickLabel',{'18Q','30Q','45Q','65Q','81Q'},...
    'YColor',[0 0 0]); %Emily Added  
xlim([1.1 2.5]);
publication_fig(0,0,1);

% figure('Name','Directional bias run length box and violin','NumberTitle','off'), hold on,
% % violin(plus_proc_times, 'xlabel',{'18Q','30Q','45Q','65Q','81Q'},'facecolor',[colour_18Q;colour_30Q;colour_45Q;colour_65Q;colour_81Q],'edgecolor','k','bw',0.1,'mc','k','medc','k--');
% violin(plus_proc_times, 'xlabel',{'30Q','45Q','65Q'},'facecolor',[colour_30Q;colour_45Q;colour_65Q],'edgecolor','k','bw',0.1,'mc','k','medc','k--');
% h2=notBoxPlot_nodotorline(plus_proc_times{1},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerSize', 5);
% h3=notBoxPlot_nodotorline(plus_proc_times{2},3,'style','line');
% set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerSize', 5);
% h4=notBoxPlot_nodotorline(plus_proc_times{3},4,'style','line');
% set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerSize', 5);
% % h5=notBoxPlot_nodotorline(plus_proc_times{4},5,'style','line');
% % set(h5.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerSize', 5);
% xlim([0.5 5.5]);
% publication_fig(0,0,1)
% axisHandle = gca; 
%     set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%     set(gca,'LineWidth',2);
%     set(gca,'FontSize',24);
%     set(gca, 'FontName', 'Arial');
%     set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
%     set(gca,'Box','on');
% ylabel('Fraction of outward processive runs');
% % ylim([-0.25 2]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off')



%% Plot CDF Plots:

figure('Name','CDF plus','NumberTitle','off'), hold on,
% h1=cdfplot(proc_plus_run_18Q);
% h1.Color=colour_18Q;
hold on, 
h2=cdfplot(proc_plus_run_30Q);
h2.Color=colour_30Q;
hold on,
h3=cdfplot(proc_plus_run_45Q);
h3.Color=colour_45Q;
hold on,
h4=cdfplot(proc_plus_run_65Q);
h4.Color=colour_65Q;
hold on,
h5=cdfplot(proc_plus_run_81Q);
h5.Color=colour_81Q;
set(gca,'xscale','log');
xlim([0 10]);
xlabel('LOG[L_{REV,+}]'); ylabel('CDF'); 
title('Plus-ended runs');
publication_fig(0,0,1);

%figure(4), hold on, 
figure('Name','CDF minus','NumberTitle','off'), hold on,
% h1=cdfplot(-1.*proc_minus_run_18Q);
% h1.Color=[0.5 0.5 0.5];
hold on, 
h2=cdfplot(-1.*proc_minus_run_30Q);
h2.Color=colour_30Q;
hold on, 
h3=cdfplot(-1.*proc_minus_run_45Q);
h3.Color=colour_45Q;
hold on,
h4=cdfplot(-1.*proc_minus_run_65Q);
h4.Color=colour_65Q;
hold on,
h5=cdfplot(-1.*proc_minus_run_81Q);
h5.Color=colour_81Q;
set(gca,'xscale','log');
xlim([0 10]);
xlabel('LOG[L_{REV,-}]'); ylabel('CDF'); 
title('Minus-ended runs');
publication_fig(0,0,1);

%% Determine Dispersion by measuring variance:
% pos_18Q=pos_all{1};
pos_30Q=pos_all{1};
pos_45Q=pos_all{2};
pos_65Q=pos_all{3};
pos_81Q=pos_all{4};

% plus_pos_18Q=pos_18Q(pos_18Q>0);
% minus_pos_18Q=pos_18Q(pos_18Q<0);

plus_pos_30Q=pos_30Q(pos_30Q>0);
minus_pos_30Q=pos_30Q(pos_30Q<0);

plus_pos_45Q=pos_45Q(pos_45Q>0);
minus_pos_45Q=pos_45Q(pos_45Q<0);

plus_pos_65Q=pos_65Q(pos_65Q>0);
minus_pos_65Q=pos_65Q(pos_65Q<0);

plus_pos_81Q=pos_81Q(pos_81Q>0);
minus_pos_81Q=pos_81Q(pos_81Q<0);

plus_pos_dist=0:0.5:11;
minus_pos_dist=-11:0.5:0;

% hist_plus_18Q=hist(plus_pos_18Q,plus_pos_dist);
hist_plus_30Q=hist(plus_pos_30Q,plus_pos_dist);
hist_plus_45Q=hist(plus_pos_45Q,plus_pos_dist);
hist_plus_65Q=hist(plus_pos_65Q,plus_pos_dist);
hist_plus_81Q=hist(plus_pos_81Q,plus_pos_dist);
% hist_minus_18Q=hist(minus_pos_18Q,minus_pos_dist);
hist_minus_30Q=hist(minus_pos_30Q,minus_pos_dist);
hist_minus_45Q=hist(minus_pos_45Q,minus_pos_dist);
hist_minus_65Q=hist(minus_pos_65Q,minus_pos_dist);
hist_minus_81Q=hist(minus_pos_81Q,minus_pos_dist);

% %figure(5), hold on, 
% figure('Name','Position frequency plus','NumberTitle','off'), hold on,
% stairs(plus_pos_dist,hist_plus_18Q./sum(hist_plus_18Q+hist_plus_30Q+hist_plus_45Q+hist_plus_65Q+hist_plus_81Q),'Color',colour_18Q,'LineWidth',2);
% hold on, 
% stairs(plus_pos_dist,hist_plus_30Q./sum(hist_plus_18Q+hist_plus_30Q+hist_plus_45Q+hist_plus_65Q+hist_plus_81Q),'Color',colour_30Q,'LineWidth',2);
% hold on, 
% stairs(plus_pos_dist,hist_plus_45Q./sum(hist_plus_18Q+hist_plus_30Q+hist_plus_45Q+hist_plus_65Q+hist_plus_81Q),'Color',colour_45Q,'LineWidth',2);
% % hold on, 
% stairs(plus_pos_dist,hist_plus_65Q./sum(hist_plus_18Q+hist_plus_30Q+hist_plus_45Q+hist_plus_65Q+hist_plus_81Q),'Color',colour_65Q,'LineWidth',2);
% hold on,
% stairs(plus_pos_dist,hist_plus_81Q./sum(hist_plus_18Q+hist_plus_30Q+hist_plus_45Q+hist_plus_65Q+hist_plus_81Q),'Color',colour_81Q,'LineWidth',2);
% set(gca, 'YScale', 'log');
% xlabel('Position (Plus)'); ylabel('Frequency'); 
% publication_fig(0,0,1);

figure('Name','Position frequency histogram','NumberTitle','off'), hold on,
% histogram(plus_pos_18Q, plus_pos_dist,'Normalization','probability','FaceColor',colour_18Q,'EdgeColor',colour_18Q,'LineWidth',2, 'FaceAlpha',0.2),hold on;
histogram(plus_pos_30Q, plus_pos_dist,'Normalization','probability','FaceColor',colour_30Q,'EdgeColor',colour_30Q,'LineWidth',2,'FaceAlpha',0.2),hold on;
histogram(plus_pos_45Q, plus_pos_dist,'Normalization','probability','FaceColor',colour_45Q,'EdgeColor',colour_45Q,'LineWidth',2,'FaceAlpha',0.2);
histogram(plus_pos_65Q, plus_pos_dist,'Normalization','probability','FaceColor',colour_65Q,'EdgeColor',colour_65Q,'LineWidth',2,'FaceAlpha',0.2);
histogram(plus_pos_81Q, plus_pos_dist,'Normalization','probability','FaceColor',colour_81Q,'EdgeColor',colour_81Q,'LineWidth',2,'FaceAlpha',0.2);
% histogram(minus_pos_18Q, minus_pos_dist,'Normalization','probability','FaceColor',colour_18Q,'EdgeColor',colour_18Q,'LineWidth',2,'FaceAlpha',0.2),hold on;
histogram(minus_pos_30Q, minus_pos_dist,'Normalization','probability','FaceColor',colour_30Q,'EdgeColor',colour_30Q,'LineWidth',2,'FaceAlpha',0.2),hold on;
histogram(minus_pos_45Q, minus_pos_dist,'Normalization','probability','FaceColor',colour_45Q,'EdgeColor',colour_45Q,'LineWidth',2,'FaceAlpha',0.2);
histogram(minus_pos_65Q, minus_pos_dist,'Normalization','probability','FaceColor',colour_65Q,'EdgeColor',colour_65Q,'LineWidth',2, 'FaceAlpha',0.2);
histogram(minus_pos_81Q, minus_pos_dist,'Normalization','probability','FaceColor',colour_81Q,'EdgeColor',colour_81Q,'LineWidth',2,'FaceAlpha',0.2);
xlabel('Position (Plus)'); ylabel('Frequency'); 
publication_fig(0,0,1);


%figure(6), hold on, 
% figure('Name','Position frequency minus','NumberTitle','off'), hold on,
% stairs(minus_pos_dist,hist_minus_18Q./sum(hist_minus_18Q+hist_minus_45Q+hist_minus_30Q+hist_minus_65Q+hist_minus_81Q),'Color',[0.5 0.5 0.5],'LineWidth',2);
% hold on, 
% stairs(minus_pos_dist,hist_minus_30Q./sum(hist_minus_18Q+hist_minus_45Q+hist_minus_30Q+hist_minus_65Q+hist_minus_81Q),'Color',colour_30Q,'LineWidth',2);
% hold on,
% stairs(minus_pos_dist,hist_minus_45Q./sum(hist_minus_18Q+hist_minus_45Q+hist_minus_30Q+hist_minus_65Q+hist_minus_81Q),'Color',colour_45Q,'LineWidth',2);
% hold on,
% stairs(minus_pos_dist,hist_minus_65Q./sum(hist_minus_18Q+hist_minus_45Q+hist_minus_30Q+hist_minus_65Q+hist_minus_81Q),'Color',colour_65Q,'LineWidth',2);
% hold on,
% stairs(minus_pos_dist,hist_minus_81Q./sum(hist_minus_18Q+hist_minus_45Q+hist_minus_30Q++hist_minus_65Q+hist_minus_81Q),'Color',colour_81Q,'LineWidth',2);
% set(gca, 'YScale', 'log');
% xlabel('Position (Minus)'); ylabel('Frequency'); 
% publication_fig(0,0,1);

%figure(1000), hold on, 
figure('Name','Dispersion','NumberTitle','off'), hold on,
% bar(1,std(plus_pos_18Q).*1000,'BarWidth',0.5,'FaceColor',colour_18Q,'Facealpha',1);
% hold on, 
% errorbar(1,std(plus_pos_18Q).*1000,0.032.*1000,'k-','LineWidth',2)
hold on,
bar(2,std(plus_pos_30Q).*1000,'BarWidth',0.5,'FaceColor',colour_30Q,'Facealpha',1);
hold on, 
errorbar(2,std(plus_pos_30Q).*1000,0.032.*1000,'k-','LineWidth',2)
hold on,
bar(3,std(plus_pos_45Q).*1000,'BarWidth',0.5,'FaceColor',colour_45Q,'Facealpha',1);
hold on, 
errorbar(3,std(plus_pos_45Q).*1000,0.059.*1000,'k-','LineWidth',2)
hold on,
bar(4,std(plus_pos_65Q).*1000,'BarWidth',0.5,'FaceColor',colour_65Q,'Facealpha',1);
hold on, 
errorbar(4,std(plus_pos_65Q).*1000,0.059.*1000,'k-','LineWidth',2)
hold on,
bar(5,std(plus_pos_81Q).*1000,'BarWidth',0.5,'FaceColor',colour_81Q,'Facealpha',1);
hold on, 
errorbar(5,std(plus_pos_81Q).*1000,0.059.*1000,'k-','LineWidth',2)
hold on,

% bar(1,-1*std(minus_pos_18Q).*1000,'BarWidth',0.5,'FaceColor',colour_18Q,'Facealpha',0.3);
% hold on, 
% errorbar(1,-1*std(minus_pos_18Q).*1000,0.057.*1000,'k-','LineWidth',2)
hold on,
bar(2,-1*std(minus_pos_30Q).*1000,'BarWidth',0.5,'FaceColor',colour_30Q,'Facealpha',0.3);
hold on, 
errorbar(2,-1*std(minus_pos_30Q).*1000,0.057.*1000,'k-','LineWidth',2)
hold on,
bar(3,-1*std(minus_pos_45Q).*1000,'BarWidth',0.5,'FaceColor',colour_45Q,'Facealpha',0.3);
hold on, 
errorbar(3,-1*std(minus_pos_45Q).*1000,0.067.*1000,'k-','LineWidth',2)
hold on,
bar(4,-1*std(minus_pos_65Q).*1000,'BarWidth',0.5,'FaceColor',colour_65Q,'Facealpha',0.3);
hold on, 
errorbar(4,-1*std(minus_pos_65Q).*1000,0.067.*1000,'k-','LineWidth',2)
hold on,
bar(5,-1*std(minus_pos_81Q).*1000,'BarWidth',0.5,'FaceColor',colour_81Q,'Facealpha',0.3);
hold on, 
errorbar(5,-1*std(minus_pos_81Q).*1000,0.067.*1000,'k-','LineWidth',2)
ylabel('Average Dispersion (nm)');
publication_fig(0,0,1);
xlim([0 5])
ylim([-7000 7000])
set(gca,'Ygrid','on');

%% Velocities:
% plus_vel_18Q=avg_vel_per_trajp{1};
plus_vel_30Q=avg_vel_per_trajp{1};
plus_vel_45Q=avg_vel_per_trajp{2};
plus_vel_65Q=avg_vel_per_trajp{3};
plus_vel_81Q=avg_vel_per_trajp{4};
% minus_vel_18Q=avg_vel_per_trajm{1};
minus_vel_30Q=avg_vel_per_trajm{1};
minus_vel_45Q=avg_vel_per_trajm{2};
minus_vel_65Q=avg_vel_per_trajm{3};
minus_vel_81Q=avg_vel_per_trajm{4};


xplus_vel=0:0.05:2;
xminus_vel=-2:0.05:0;

% xhist_plus_vel_18Q=hist(avg_vel_per_trajp{1},xplus_vel);
% xhist_minus_vel_18Q=hist(avg_vel_per_trajm{1},xminus_vel);
xhist_plus_vel_30Q=hist(avg_vel_per_trajp{1},xplus_vel);
xhist_minus_vel_30Q=hist(avg_vel_per_trajm{1},xminus_vel);
xhist_plus_vel_45Q=hist(avg_vel_per_trajp{2},xplus_vel);
xhist_minus_vel_45Q=hist(avg_vel_per_trajm{2},xminus_vel);
xhist_plus_vel_65Q=hist(avg_vel_per_trajp{3},xplus_vel);
xhist_minus_vel_65Q=hist(avg_vel_per_trajm{3},xminus_vel);
xhist_plus_vel_81Q=hist(avg_vel_per_trajp{4},xplus_vel);
xhist_minus_vel_81Q=hist(avg_vel_per_trajm{4},xminus_vel);

plus_vel_nonan_30Q=avg_vel_per_trajp{1}(find(~isnan(avg_vel_per_trajp{1})));
plus_vel_nonan_45Q=avg_vel_per_trajp{2}(find(~isnan(avg_vel_per_trajp{2})));
plus_vel_nonan_65Q=avg_vel_per_trajp{3}(find(~isnan(avg_vel_per_trajp{3})));
plus_vel_nonan_81Q=avg_vel_per_trajp{4}(find(~isnan(avg_vel_per_trajp{4})));

minus_vel_nonan_30Q=avg_vel_per_trajm{1}(find(~isnan(avg_vel_per_trajm{1})));
minus_vel_nonan_45Q=avg_vel_per_trajm{2}(find(~isnan(avg_vel_per_trajm{2})));
minus_vel_nonan_65Q=avg_vel_per_trajm{3}(find(~isnan(avg_vel_per_trajm{3})));
minus_vel_nonan_81Q=avg_vel_per_trajm{4}(find(~isnan(avg_vel_per_trajm{4})));

plus_lr_nonan_30Q=avg_lr_per_trajp{1}(find(~isnan(avg_lr_per_trajp{1})));
plus_lr_nonan_45Q=avg_lr_per_trajp{2}(find(~isnan(avg_lr_per_trajp{2})));
plus_lr_nonan_65Q=avg_lr_per_trajp{3}(find(~isnan(avg_lr_per_trajp{3})));
plus_lr_nonan_81Q=avg_lr_per_trajp{4}(find(~isnan(avg_lr_per_trajp{4})));

minus_lr_nonan_30Q=avg_lr_per_trajm{1}(find(~isnan(avg_lr_per_trajm{1})));
minus_lr_nonan_45Q=avg_lr_per_trajm{2}(find(~isnan(avg_lr_per_trajm{2})));
minus_lr_nonan_65Q=avg_lr_per_trajm{3}(find(~isnan(avg_lr_per_trajm{3})));
minus_lr_nonan_81Q=avg_lr_per_trajm{4}(find(~isnan(avg_lr_per_trajm{4})));

SEM_plus_vel_30Q=std(plus_vel_nonan_30Q)/sqrt(length(plus_vel_nonan_30Q));
SEM_plus_vel_45Q=std(plus_vel_nonan_45Q)/sqrt(length(plus_vel_nonan_45Q));
SEM_plus_vel_65Q=std(plus_vel_nonan_65Q)/sqrt(length(plus_vel_nonan_65Q));
SEM_plus_vel_81Q=std(plus_vel_nonan_81Q)/sqrt(length(plus_vel_nonan_81Q));

SEM_minus_vel_30Q=std(minus_vel_nonan_30Q)/sqrt(length(minus_vel_nonan_30Q));
SEM_minus_vel_45Q=std(minus_vel_nonan_45Q)/sqrt(length(minus_vel_nonan_45Q));
SEM_minus_vel_65Q=std(minus_vel_nonan_65Q)/sqrt(length(minus_vel_nonan_65Q));
SEM_minus_vel_81Q=std(minus_vel_nonan_81Q)/sqrt(length(minus_vel_nonan_81Q));

SEM_plus_lr_30Q=std(plus_lr_nonan_30Q)/sqrt(length(plus_lr_nonan_30Q));
SEM_plus_lr_45Q=std(plus_lr_nonan_45Q)/sqrt(length(plus_lr_nonan_45Q));
SEM_plus_lr_65Q=std(plus_lr_nonan_65Q)/sqrt(length(plus_lr_nonan_65Q));
SEM_plus_lr_81Q=std(plus_lr_nonan_81Q)/sqrt(length(plus_lr_nonan_81Q));

SEM_minus_lr_30Q=std(minus_lr_nonan_30Q)/sqrt(length(minus_lr_nonan_30Q));
SEM_minus_lr_45Q=std(minus_lr_nonan_45Q)/sqrt(length(minus_lr_nonan_45Q));
SEM_minus_lr_65Q=std(minus_lr_nonan_65Q)/sqrt(length(minus_lr_nonan_65Q));
SEM_minus_lr_81Q=std(minus_lr_nonan_81Q)/sqrt(length(minus_lr_nonan_81Q));

%% Travel Distance vs Velocity (um/s)

bin_interval_vel_p=0:0.1:2.5;
bin_interval_vel_m=-2.5:0.1:0;

figure('Name','Processive Velocity distribution','NumberTitle','off'), hold on,
histogram(avg_vel_per_trajp{1},bin_interval_vel_p,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajp{2},bin_interval_vel_p,'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajp{3},bin_interval_vel_p,'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajp{4},bin_interval_vel_p,'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajm{1},bin_interval_vel_m,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajm{2},bin_interval_vel_m,'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajm{3},bin_interval_vel_m,'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_vel_per_trajm{4},bin_interval_vel_m,'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
xlabel('Velocity (\mum/s)');
ylabel('Frequency');
publication_fig(0,0,1);


bin_interval=0:0.25:10;
bin_interval_minus=-10:0.25:0;

figure('Name','Run length distribution','NumberTitle','off'), hold on,
histogram(avg_lr_per_trajp{1},bin_interval,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajp{2},bin_interval,'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajp{3},bin_interval,'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajp{4},bin_interval,'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajm{1},bin_interval_minus,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajm{2},bin_interval_minus,'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajm{3},bin_interval_minus,'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(avg_lr_per_trajm{4},bin_interval_minus,'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
xlabel('Run length (\mum)');
ylabel('Frequency');
publication_fig(0,0,1);

% all_err_bars=[SEM_plus_vel_30Q SEM_plus_vel_45Q SEM_plus_vel_65Q SEM_plus_vel_81Q; SEM_minus_vel_30Q SEM_minus_vel_45Q SEM_minus_vel_65Q SEM_minus_vel_81Q];
% values=[avg_vel_per_trajp{1} avg_vel_per_trajp{2} avg_vel_per_trajp{3} avg_vel_per_trajm{4}; avg_vel_per_trajm{1} avg_vel_per_trajm{2} avg_vel_per_trajm{3} avg_vel_per_trajm{4}];
% figure('Name','Velocities Bar','NumberTitle','off'), hold on,
% b2=bar(values,'grouped');
% b2(1).FaceColor=colour_30Q;
% b2(2).FaceColor=colour_45Q;
% b2(3).FaceColor=colour_65Q;
% b2(4).FaceColor=colour_81Q;
% b2(1).FaceAlpha=0.5;
% b2(2).FaceAlpha=0.5;
% b2(3).FaceAlpha=0.5;
% b2(4).FaceAlpha=0.5;
% hold on
% [ngroups,nbars] = size(values);
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, values(:,i), all_err_bars(:,i), 'k', 'linestyle', 'none','LineWidth',2);
% end
% xticks([1 2]);
% xticklabels({'Anterograde', 'Retrograde'});
% ylabel('Mean Velocity \mum/s');
% publication_fig(0,0,1)


%thresholding the data to plot a narrower range in y, and plot all values
%above a certain threshold at one spot on the y axis for > threshold.
%Untested 20230923.

vel_threshold=2.5;
lr_threshold=25;

plus_vel_thres{1}=plus_vel_nonan_30Q(find(plus_vel_nonan_30Q<=vel_threshold));
plus_vel_thres{2}=plus_vel_nonan_45Q(plus_vel_nonan_45Q<=vel_threshold);
plus_vel_thres{3}=plus_vel_nonan_65Q(find(plus_vel_nonan_65Q<=vel_threshold));
plus_vel_thres{4}=plus_vel_nonan_81Q(find(plus_vel_nonan_81Q<=vel_threshold));

minus_vel_thres{1}=minus_vel_nonan_30Q(find(-1*minus_vel_nonan_30Q<=vel_threshold));
minus_vel_thres{2}=minus_vel_nonan_45Q(find(-1*minus_vel_nonan_45Q<=vel_threshold));
minus_vel_thres{3}=minus_vel_nonan_65Q(find(-1*minus_vel_nonan_65Q<=vel_threshold));
minus_vel_thres{4}=minus_vel_nonan_81Q(find(-1*minus_vel_nonan_81Q<=vel_threshold));

plus_lr_thres{1}=plus_lr_nonan_30Q(find(plus_lr_nonan_30Q<=lr_threshold));
plus_lr_thres{2}=plus_lr_nonan_45Q(find(plus_lr_nonan_45Q<=lr_threshold));
plus_lr_thres{3}=plus_lr_nonan_65Q(find(plus_lr_nonan_65Q<=lr_threshold));
plus_lr_thres{4}=plus_lr_nonan_81Q(find(plus_lr_nonan_81Q<=lr_threshold));

minus_lr_thres{1}=minus_lr_nonan_30Q(find(-1*minus_lr_nonan_30Q<=lr_threshold));
minus_lr_thres{2}=minus_lr_nonan_45Q(find(-1*minus_lr_nonan_45Q<=lr_threshold));
minus_lr_thres{3}=minus_lr_nonan_65Q(find(-1*minus_lr_nonan_65Q<=lr_threshold));
minus_lr_thres{4}=minus_lr_nonan_81Q(find(-1*minus_lr_nonan_81Q<=lr_threshold));

plus_vel_above_thres{1}=plus_vel_nonan_30Q(find(plus_vel_nonan_30Q>vel_threshold));
plus_vel_above_thres{2}=plus_vel_nonan_45Q(find(plus_vel_nonan_45Q>vel_threshold));
plus_vel_above_thres{3}=plus_vel_nonan_65Q(find(plus_vel_nonan_65Q>vel_threshold));
plus_vel_above_thres{4}=plus_vel_nonan_81Q(find(plus_vel_nonan_81Q>vel_threshold));

minus_vel_above_thres{1}=minus_vel_nonan_30Q(find(-1*minus_vel_nonan_30Q>vel_threshold));
minus_vel_above_thres{2}=minus_vel_nonan_45Q(find(-1*minus_vel_nonan_45Q>vel_threshold));
minus_vel_above_thres{3}=minus_vel_nonan_65Q(find(-1*minus_vel_nonan_65Q>vel_threshold));
minus_vel_above_thres{4}=minus_vel_nonan_81Q(find(-1*minus_vel_nonan_81Q>vel_threshold));

plus_lr_above_thres{1}=plus_lr_nonan_30Q(find(plus_lr_nonan_30Q>lr_threshold));
plus_lr_above_thres{2}=plus_lr_nonan_45Q(find(plus_lr_nonan_45Q>lr_threshold));
plus_lr_above_thres{3}=plus_lr_nonan_65Q(find(plus_lr_nonan_65Q>lr_threshold));
plus_lr_above_thres{4}=plus_lr_nonan_81Q(find(plus_lr_nonan_81Q>lr_threshold));

minus_lr_above_thres{1}=minus_lr_nonan_30Q(find(-1*minus_lr_nonan_30Q>lr_threshold));
minus_lr_above_thres{2}=minus_lr_nonan_45Q(find(-1*minus_lr_nonan_45Q>lr_threshold));
minus_lr_above_thres{3}=minus_lr_nonan_65Q(find(-1*minus_lr_nonan_65Q>lr_threshold));
minus_lr_above_thres{4}=minus_lr_nonan_81Q(find(-1*minus_lr_nonan_81Q>lr_threshold));

n_p_vel_above_thres{1}=numel(plus_vel_above_thres{1});
n_p_vel_above_thres{2}=numel(plus_vel_above_thres{2});
n_p_vel_above_thres{3}=numel(plus_vel_above_thres{3});
n_p_vel_above_thres{4}=numel(plus_vel_above_thres{4});

n_m_vel_above_thres{1}=numel(minus_vel_above_thres{1});
n_m_vel_above_thres{2}=numel(minus_vel_above_thres{2});
n_m_vel_above_thres{3}=numel(minus_vel_above_thres{3});
n_m_vel_above_thres{4}=numel(minus_vel_above_thres{4});

n_p_lr_above_thres{1}=numel(plus_lr_above_thres{1});
n_p_lr_above_thres{2}=numel(plus_lr_above_thres{2});
n_p_lr_above_thres{3}=numel(plus_lr_above_thres{3});
n_p_lr_above_thres{4}=numel(plus_lr_above_thres{4});

n_m_lr_above_thres{1}=numel(minus_lr_above_thres{1});
n_m_lr_above_thres{2}=numel(minus_lr_above_thres{2});
n_m_lr_above_thres{3}=numel(minus_lr_above_thres{3});
n_m_lr_above_thres{4}=numel(minus_lr_above_thres{4});

figure('Name','Thresholded Velocities','NumberTitle','off'), hold on,
plot(1,mean(plus_vel_nonan_30Q),'k_','MarkerSize',30);
plot(2,mean(plus_vel_nonan_45Q),'k_','MarkerSize',30);
plot(3,mean(plus_vel_nonan_65Q),'k_','MarkerSize',30);
plot(4,mean(plus_vel_nonan_81Q),'k_','MarkerSize',30);
plot(6,-1*mean(minus_vel_nonan_30Q),'k_','MarkerSize',30);
plot(7,-1*mean(minus_vel_nonan_45Q),'k_','MarkerSize',30);
plot(8,-1*mean(minus_vel_nonan_65Q),'k_','MarkerSize',30);
plot(9,-1*mean(minus_vel_nonan_81Q),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(plus_vel_thres{1}),1),plus_vel_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_vel_thres{2}),1),plus_vel_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_vel_thres{3}),1),plus_vel_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_vel_thres{4}),1),plus_vel_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_vel_thres{1}),1),-1*minus_vel_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_vel_thres{2}),1),-1*minus_vel_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_vel_thres{3}),1),-1*minus_vel_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_vel_thres{4}),1),-1*minus_vel_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(plus_vel_above_thres{1}),1),(vel_threshold+0.2)*ones(n_p_vel_above_thres{1},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_vel_above_thres{2}),1),(vel_threshold+0.2)*ones(n_p_vel_above_thres{2},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_vel_above_thres{3}),1),(vel_threshold+0.2)*ones(n_p_vel_above_thres{3},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_vel_above_thres{4}),1),(vel_threshold+0.2)*ones(n_p_vel_above_thres{4},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_vel_above_thres{1}),1),(vel_threshold+0.2)*ones(n_m_vel_above_thres{1},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_vel_above_thres{2}),1),(vel_threshold+0.2)*ones(n_m_vel_above_thres{2},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_vel_above_thres{3}),1),(vel_threshold+0.2)*ones(n_m_vel_above_thres{3},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_vel_above_thres{4}),1),(vel_threshold+0.2)*ones(n_m_vel_above_thres{4},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(plus_vel_nonan_30Q),SEM_plus_vel_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(plus_vel_nonan_45Q),SEM_plus_vel_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(plus_vel_nonan_65Q),SEM_plus_vel_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(plus_vel_nonan_81Q),SEM_plus_vel_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
errorbar(6,mean(-1*minus_vel_nonan_30Q),-1*SEM_minus_vel_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(7,mean(-1*minus_vel_nonan_45Q),-1*SEM_minus_vel_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(8,mean(-1*minus_vel_nonan_65Q),-1*SEM_minus_vel_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(9,mean(-1*minus_vel_nonan_81Q),-1*SEM_minus_vel_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 10]);
ylabel('Velocity (\mum/s)');
xticks([2.5 7.5]);
xticklabels({'Ant', 'Ret'});
publication_fig(0,0,1);
box('off');
set(gca, 'FontSize', 40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Thresholded Run Lengths','NumberTitle','off'), hold on,
plot(1,mean(plus_lr_nonan_30Q),'k_','MarkerSize',30);
plot(2,mean(plus_lr_nonan_45Q),'k_','MarkerSize',30);
plot(3,mean(plus_lr_nonan_65Q),'k_','MarkerSize',30);
plot(4,mean(plus_lr_nonan_81Q),'k_','MarkerSize',30);
plot(6,-1*mean(minus_lr_nonan_30Q),'k_','MarkerSize',30);
plot(7,-1*mean(minus_lr_nonan_45Q),'k_','MarkerSize',30);
plot(8,-1*mean(minus_lr_nonan_65Q),'k_','MarkerSize',30);
plot(9,-1*mean(minus_lr_nonan_81Q),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(plus_lr_thres{1}),1),plus_lr_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_lr_thres{2}),1),plus_lr_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_lr_thres{3}),1),plus_lr_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_lr_thres{4}),1),plus_lr_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_lr_thres{1}),1),-1*minus_lr_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_lr_thres{2}),1),-1*minus_lr_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_lr_thres{3}),1),-1*minus_lr_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_lr_thres{4}),1),-1*minus_lr_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(plus_lr_above_thres{1}),1),(lr_threshold+1)*ones(n_p_lr_above_thres{1},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_lr_above_thres{2}),1),(lr_threshold+1)*ones(n_p_lr_above_thres{2},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_lr_above_thres{3}),1),(lr_threshold+1)*ones(n_p_lr_above_thres{3},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_lr_above_thres{4}),1),(lr_threshold+1)*ones(n_p_lr_above_thres{4},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_lr_above_thres{1}),1),(lr_threshold+1)*ones(n_m_lr_above_thres{1},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_lr_above_thres{2}),1),(lr_threshold+1)*ones(n_m_lr_above_thres{2},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_lr_above_thres{3}),1),(lr_threshold+1)*ones(n_m_lr_above_thres{3},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_lr_above_thres{4}),1),(lr_threshold+1)*ones(n_m_lr_above_thres{4},1),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(plus_lr_nonan_30Q),SEM_plus_lr_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(plus_lr_nonan_45Q),SEM_plus_lr_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(plus_lr_nonan_65Q),SEM_plus_lr_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(plus_lr_nonan_81Q),SEM_plus_lr_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
errorbar(6,mean(-1*minus_lr_nonan_30Q),-1*SEM_minus_lr_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(7,mean(-1*minus_lr_nonan_45Q),-1*SEM_minus_lr_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(8,mean(-1*minus_lr_nonan_65Q),-1*SEM_minus_lr_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(9,mean(-1*minus_lr_nonan_81Q),-1*SEM_minus_lr_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 10]);
xticks([2.5 7.5]);
xticklabels({'Ant', 'Ret'});
ylabel('Run Length (\mum)');
publication_fig(0,0,1);
box('off');
set(gca, 'FontSize', 40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');


figure('Name','Velocities','NumberTitle','off'), hold on,
plot(1,mean(plus_vel_nonan_30Q),'k_','MarkerSize',30);
plot(2,mean(plus_vel_nonan_45Q),'k_','MarkerSize',30);
plot(3,mean(plus_vel_nonan_65Q),'k_','MarkerSize',30);
plot(4,mean(plus_vel_nonan_81Q),'k_','MarkerSize',30);
plot(6,-1*mean(minus_vel_nonan_30Q),'k_','MarkerSize',30);
plot(7,-1*mean(minus_vel_nonan_45Q),'k_','MarkerSize',30);
plot(8,-1*mean(minus_vel_nonan_65Q),'k_','MarkerSize',30);
plot(9,-1*mean(minus_vel_nonan_81Q),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(plus_vel_nonan_30Q),1),plus_vel_nonan_30Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_vel_nonan_45Q),1),plus_vel_nonan_45Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_vel_nonan_65Q),1),plus_vel_nonan_65Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_vel_nonan_81Q),1),plus_vel_nonan_81Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_vel_nonan_30Q),1),-1*minus_vel_nonan_30Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_vel_nonan_45Q),1),-1*minus_vel_nonan_45Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_vel_nonan_65Q),1),-1*minus_vel_nonan_65Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_vel_nonan_81Q),1),-1*minus_vel_nonan_81Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(plus_vel_nonan_30Q),SEM_plus_vel_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(plus_vel_nonan_45Q),SEM_plus_vel_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(plus_vel_nonan_65Q),SEM_plus_vel_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(plus_vel_nonan_81Q),SEM_plus_vel_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
errorbar(6,mean(-1*minus_vel_nonan_30Q),-1*SEM_minus_vel_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(7,mean(-1*minus_vel_nonan_45Q),-1*SEM_minus_vel_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(8,mean(-1*minus_vel_nonan_65Q),-1*SEM_minus_vel_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(9,mean(-1*minus_vel_nonan_81Q),-1*SEM_minus_vel_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 10]);
ylabel('Velocity (\mum/s)');
publication_fig(0,0,1);
box('off');
set(gca, 'FontSize', 40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Run Lengths','NumberTitle','off'), hold on,
plot(1,mean(plus_lr_nonan_30Q),'k_','MarkerSize',30);
plot(2,mean(plus_lr_nonan_45Q),'k_','MarkerSize',30);
plot(3,mean(plus_lr_nonan_65Q),'k_','MarkerSize',30);
plot(4,mean(plus_lr_nonan_81Q),'k_','MarkerSize',30);
plot(6,-1*mean(minus_lr_nonan_30Q),'k_','MarkerSize',30);
plot(7,-1*mean(minus_lr_nonan_45Q),'k_','MarkerSize',30);
plot(8,-1*mean(minus_lr_nonan_65Q),'k_','MarkerSize',30);
plot(9,-1*mean(minus_lr_nonan_81Q),'k_','MarkerSize',30);
errorbar(1,mean(plus_lr_nonan_30Q),SEM_plus_lr_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(plus_lr_nonan_45Q),SEM_plus_lr_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(plus_lr_nonan_65Q),SEM_plus_lr_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(plus_lr_nonan_81Q),SEM_plus_lr_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
errorbar(6,mean(-1*minus_lr_nonan_30Q),-1*SEM_minus_lr_30Q,'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(7,mean(-1*minus_lr_nonan_45Q),-1*SEM_minus_lr_45Q,'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(8,mean(-1*minus_lr_nonan_65Q),-1*SEM_minus_lr_65Q,'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(9,mean(-1*minus_lr_nonan_81Q),-1*SEM_minus_lr_81Q,'Color',colour_81Q,'LineWidth',4,'CapSize',0);
scatter(1-0.125+0.25*rand(numel(plus_lr_nonan_30Q),1),plus_lr_nonan_30Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(plus_lr_nonan_45Q),1),plus_lr_nonan_45Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(plus_lr_nonan_65Q),1),plus_lr_nonan_65Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(plus_lr_nonan_81Q),1),plus_lr_nonan_81Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(minus_lr_nonan_30Q),1),-1*minus_lr_nonan_30Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(minus_lr_nonan_45Q),1),-1*minus_lr_nonan_45Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(minus_lr_nonan_65Q),1),-1*minus_lr_nonan_65Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(minus_lr_nonan_81Q),1),-1*minus_lr_nonan_81Q,50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
xlim([0 10]);
ylabel('Run Length (\mum)');
publication_fig(0,0,1);
box('off');
set(gca, 'FontSize', 40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');


%figure(8), hold on, subplot(1,2,1)
figure('Name','Run length velocity plus and minus','NumberTitle','off'), hold on, subplot(1,2,1)
scatter(avg_lr_per_trajp{1},avg_vel_per_trajp{1},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(avg_lr_per_trajp{2},avg_vel_per_trajp{2},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(avg_lr_per_trajp{3},avg_vel_per_trajp{3},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(avg_lr_per_trajp{4},avg_vel_per_trajp{4},'o','MarkerFaceColor',colour_81Q,'MarkerEdgeColor',colour_81Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 25]); ylim([0 3.5]);
    
hold on, 
subplot(1,2,2)
scatter(-1*avg_lr_per_trajm{1},-1*avg_vel_per_trajm{1},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(-1*avg_lr_per_trajm{2},-1*avg_vel_per_trajm{2},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(-1*avg_lr_per_trajm{3},-1*avg_vel_per_trajm{3},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
scatter(-1*avg_lr_per_trajm{4},-1*avg_vel_per_trajm{4},'o','MarkerFaceColor',colour_81Q,'MarkerEdgeColor',colour_81Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2),hold on;
xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 25]); ylim([0 3.5]);

%testing difference heatmaps of 30Q vs 81Q

rl_81Q_no_nan=avg_lr_per_trajp{4}(find(~isnan(avg_lr_per_trajp{4})));
rl_30Q_no_nan=avg_lr_per_trajp{1}(find(~isnan(avg_lr_per_trajp{1})));
vel_81Q_no_nan=avg_vel_per_trajp{4}(find(~isnan(avg_vel_per_trajp{4})));
vel_30Q_no_nan=avg_vel_per_trajp{1}(find(~isnan(avg_vel_per_trajp{1})));

samp_81Q_rl=randsample(rl_81Q_no_nan, 100, true); %taking a sample of 100 with replacement to get an equal number for both conditions
samp_30Q_rl=randsample(rl_30Q_no_nan, 100, true); %taking a sample of 100 with replacement to get an equal number for both conditions
samp_81Q_vel=randsample(vel_81Q_no_nan, 100, true); %taking a sample of 100 with replacement to get an equal number for both conditions
samp_30Q_vel=randsample(vel_30Q_no_nan, 100, true); %taking a sample of 100 with replacement to get an equal number for both conditions

Nbins = 25;
x=samp_81Q_rl-samp_30Q_rl;
y=samp_81Q_vel-samp_30Q_vel;

figure, plot(x,y,'ko')
hold on, hb = binscatter(x,y);
hb.NumBins = [Nbins Nbins]; %[x, y]
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.5;
ax = gca;
colormap(ax,'parula')
xlabel('Run length Difference 81Q-30Q (\mum)');
ylabel('Velocity Difference 81Q-30Q(\mum/s)');

% Calculate the shift distances
shift_distances = sqrt((samp_81Q_rl-samp_30Q_rl).^2 + (samp_81Q_vel-samp_30Q_vel).^2);
% Determine the direction of the shift (towards or away from the origin)
original_distances = sqrt(samp_30Q_rl.^2 + samp_30Q_vel.^2);
new_distances = sqrt(samp_81Q_rl.^2 + samp_81Q_vel.^2);
shift_direction = new_distances - original_distances;
% Create a colormap based on the shift distances and direction
shift_intensity = min(max(shift_distances / max(shift_distances), 0), 1); % Normalize
colors = zeros(length(shift_intensity), 3);
colors(shift_direction > 0, 1) = shift_intensity(shift_direction > 0); % Red for away
colors(shift_direction <= 0, 2) = shift_intensity(shift_direction <= 0); % Green for towards
% Plot the original points
figure,
scatter(samp_30Q_rl, samp_30Q_vel, 36, 'k', 'filled');
hold on
% Overlay the shifted points with the colormap
scatter(samp_81Q_rl, samp_81Q_vel, 36, colors, 'filled');
% Add a colorbar if needed
colormap([linspace(0, 1, 256)', linspace(1, 0, 256)', zeros(256, 1)]); % Red to green colormap
c = colorbar;
c.Label.String = 'Shift Intensity';
xlabel('Run length Shift 81Q-30Q (\mum)');
ylabel('Velocity Shift 81Q-30Q(\mum/s)');


%plotting the distributions of the products from run length and velocity
prod_rl_vel_30Q=avg_vel_per_trajp{1}.*avg_lr_per_trajp{1};
prod_rl_vel_45Q=avg_vel_per_trajp{2}.*avg_lr_per_trajp{2};
prod_rl_vel_65Q=avg_vel_per_trajp{3}.*avg_lr_per_trajp{3};
prod_rl_vel_81Q=avg_vel_per_trajp{4}.*avg_lr_per_trajp{4};

prod_bins=0:0.5:20;
% figure('Name','Run length Velocity Product Histogram','NumberTitle','off'), hold on
figure(1001), hold on
histogram(prod_rl_vel_30Q, prod_bins, 'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_45Q, prod_bins, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_65Q, prod_bins, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_81Q, prod_bins, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');

options = statset('Display','final','MaxIter',1000,'TolX',1e-8,'robust','on'); %options for the gaussian mixture distribution fit

for kdist=1 %The number of gaussian fits to perform, this can be changed to a range and optimized by minimizing the BIC criteria
fit_prod_30Q = gmdistribution.fit(prod_rl_vel_30Q,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_45Q = gmdistribution.fit(prod_rl_vel_45Q,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_65Q = gmdistribution.fit(prod_rl_vel_65Q,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_81Q = gmdistribution.fit(prod_rl_vel_81Q,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions

bic_prod_30Q(kdist)=fit_prod_30Q.BIC;
bic_prod_45Q(kdist)=fit_prod_45Q.BIC;
bic_prod_65Q(kdist)=fit_prod_65Q.BIC;
bic_prod_81Q(kdist)=fit_prod_81Q.BIC;

pdf_30Q_fit=pdf(fit_prod_30Q,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_45Q_fit=pdf(fit_prod_45Q,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_65Q_fit=pdf(fit_prod_65Q,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_81Q_fit=pdf(fit_prod_81Q,prod_bins');  %Probability density function for the gaussian fit (for plotting)

pdf_30Q_fit_norm=pdf_30Q_fit./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_45Q_fit_norm=pdf_45Q_fit./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_65Q_fit_norm=pdf_65Q_fit./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_81Q_fit_norm=pdf_81Q_fit./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)

    
figure(10)
plot(bic_prod_30Q,'b--');
plot(bic_prod_45Q,'b--');
plot(bic_prod_65Q,'b--');
plot(bic_prod_81Q,'b--');
    
figure(1000), hold on
plot(prod_bins,0.035*pdf_30Q_fit,'Color',colour_30Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_45Q_fit,'Color',colour_45Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_65Q_fit,'Color',colour_65Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_81Q_fit,'Color',colour_81Q,'LineStyle','-','LineWidth', 3);
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');

figure(1001), hold on
plot(prod_bins,pdf_30Q_fit_norm,'Color',colour_30Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_45Q_fit_norm,'Color',colour_45Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_65Q_fit_norm,'Color',colour_65Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_81Q_fit_norm,'Color',colour_81Q,'LineStyle','-','LineWidth', 3);
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');
end

prod_rl_vel_30Q_m=avg_vel_per_trajm{1}.*avg_lr_per_trajm{1};
prod_rl_vel_45Q_m=avg_vel_per_trajm{2}.*avg_lr_per_trajm{2};
prod_rl_vel_65Q_m=avg_vel_per_trajm{3}.*avg_lr_per_trajm{3};
prod_rl_vel_81Q_m=avg_vel_per_trajm{4}.*avg_lr_per_trajm{4};

% figure('Name','Run length Velocity Product Histogram Minus','NumberTitle','off'), hold on
figure(2001), hold on
histogram(prod_rl_vel_30Q_m, prod_bins, 'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_45Q_m, prod_bins, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_65Q_m, prod_bins, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
histogram(prod_rl_vel_81Q_m, prod_bins, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.5,'Normalization','probability');
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');

options = statset('Display','final','MaxIter',1000,'TolX',1e-8,'robust','on'); %options for the gaussian mixture distribution fit

for kdist=1 %The number of gaussian fits to perform, this can be changed to a range and optimized by minimizing the BIC criteria
fit_prod_30Q_m = gmdistribution.fit(prod_rl_vel_30Q_m,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_45Q_m = gmdistribution.fit(prod_rl_vel_45Q_m,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_65Q_m = gmdistribution.fit(prod_rl_vel_65Q_m,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
fit_prod_81Q_m = gmdistribution.fit(prod_rl_vel_81Q_m,kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions

bic_prod_30Q_m(kdist)=fit_prod_30Q_m.BIC;
bic_prod_45Q_m(kdist)=fit_prod_45Q_m.BIC;
bic_prod_65Q_m(kdist)=fit_prod_65Q_m.BIC;
bic_prod_81Q_m(kdist)=fit_prod_81Q_m.BIC;

pdf_30Q_fit_m=pdf(fit_prod_30Q_m,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_45Q_fit_m=pdf(fit_prod_45Q_m,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_65Q_fit_m=pdf(fit_prod_65Q_m,prod_bins');  %Probability density function for the gaussian fit (for plotting)
pdf_81Q_fit_m=pdf(fit_prod_81Q_m,prod_bins');  %Probability density function for the gaussian fit (for plotting)

pdf_30Q_fit_norm_m=pdf_30Q_fit_m./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_45Q_fit_norm_m=pdf_45Q_fit_m./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_65Q_fit_norm_m=pdf_65Q_fit_m./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
pdf_81Q_fit_norm_m=pdf_81Q_fit_m./(kdist*numel(prod_bins));  %Probability density function for the gaussian fit (for plotting)
    
figure(20)
plot(bic_prod_30Q_m,'b--');
plot(bic_prod_45Q_m,'b--');
plot(bic_prod_65Q_m,'b--');
plot(bic_prod_81Q_m,'b--');
    
figure(2000), hold on
plot(prod_bins,0.035*pdf_30Q_fit_m,'Color',colour_30Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_45Q_fit_m,'Color',colour_45Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_65Q_fit_m,'Color',colour_65Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,0.035*pdf_81Q_fit_m,'Color',colour_81Q,'LineStyle','-','LineWidth', 3);
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');

figure(2001), hold on
plot(prod_bins,pdf_30Q_fit_norm_m,'Color',colour_30Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_45Q_fit_norm_m,'Color',colour_45Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_65Q_fit_norm_m,'Color',colour_65Q,'LineStyle','-','LineWidth', 3);
plot(prod_bins,pdf_81Q_fit_norm_m,'Color',colour_81Q,'LineStyle','-','LineWidth', 3);
xlabel('Run length*Velocity(\mum^2/s)');
ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','off');
end 
% figure('Name','Run length velocity plus','NumberTitle','off'), hold on
% scatter(avg_lr_per_trajp{1},avg_vel_per_trajp{1},'o','MarkerFaceColor',colour_18Q,'MarkerEdgeColor',colour_18Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(avg_lr_per_trajp{2},avg_vel_per_trajp{2},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(avg_lr_per_trajp{3},avg_vel_per_trajp{3},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(avg_lr_per_trajp{4},avg_vel_per_trajp{4},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% set(gca,'LineWidth',2);
% set(gca,'FontSize',24);
% set(gca, 'FontName', 'Arial');
% set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
% set(gca,'Box','on');
% xlim([0 25]); ylim([0 10]);
% 
% figure('Name','Run length velocity minus','NumberTitle','off'), hold on
% scatter(-1*avg_lr_per_trajm{1},-1*avg_vel_per_trajm{1},'o','MarkerFaceColor',colour_18Q,'MarkerEdgeColor',colour_18Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(-1*avg_lr_per_trajm{2},-1*avg_vel_per_trajm{2},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(-1*avg_lr_per_trajm{3},-1*avg_vel_per_trajm{3},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(-1*avg_lr_per_trajm{4},-1*avg_vel_per_trajm{4},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% set(gca,'LineWidth',2);
% set(gca,'FontSize',24);
% set(gca, 'FontName', 'Arial');
% set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
% set(gca,'Box','on');
% xlim([0 25]); ylim([0 10]);
    

figure('Name','Run length velocity','NumberTitle','off'), hold on,
% scatter(avg_lr_per_trajp{1},avg_vel_per_trajp{1},'o','MarkerFaceColor',colour_18Q,'MarkerEdgeColor',colour_18Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(avg_lr_per_trajp{1},avg_vel_per_trajp{1},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(avg_lr_per_trajp{2},avg_vel_per_trajp{2},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(avg_lr_per_trajp{3},avg_vel_per_trajp{3},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(avg_lr_per_trajp{4},avg_vel_per_trajp{4},'o','MarkerFaceColor',colour_81Q,'MarkerEdgeColor',colour_81Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% scatter(-1*avg_lr_per_trajm{1},avg_vel_per_trajm{1},'o','MarkerFaceColor',colour_18Q,'MarkerEdgeColor',colour_18Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(-1*avg_lr_per_trajm{1},avg_vel_per_trajm{1},'o','MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(-1*avg_lr_per_trajm{2},avg_vel_per_trajm{2},'o','MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(-1*avg_lr_per_trajm{3},avg_vel_per_trajm{3},'o','MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(-1*avg_lr_per_trajm{4},avg_vel_per_trajm{4},'o','MarkerFaceColor',colour_81Q,'MarkerEdgeColor',colour_81Q,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 25]); ylim([-7 7]);

% rev_rate_18Q=rev_rates{1}(rev_rates{1}>0);
rev_rate_30Q=rev_rates{1}(rev_rates{1}>0);
rev_rate_45Q=rev_rates{2}(rev_rates{2}>0);
rev_rate_65Q=rev_rates{3}(rev_rates{3}>0);
rev_rate_81Q=rev_rates{4}(rev_rates{4}>0);

%Reversal rates
figure('Name','Reversal rates','NumberTitle','off'), hold on,
% histogram(rev_rate_18Q,0:0.01:1,'Normalization','probability','FaceColor',colour_18Q,'EdgeColor',colour_18Q,'LineWidth',2,'FaceAlpha',.2);
histogram(rev_rate_30Q,0:0.01:1,'Normalization','probability','FaceColor',colour_30Q,'EdgeColor',colour_30Q,'LineWidth',2,'FaceAlpha',.2);
histogram(rev_rate_45Q,0:0.01:1,'Normalization','probability','FaceColor',colour_45Q,'EdgeColor',colour_45Q,'LineWidth',2,'FaceAlpha',.2);
histogram(rev_rate_65Q,0:0.01:1,'Normalization','probability','FaceColor',colour_65Q,'EdgeColor',colour_65Q,'LineWidth',2,'FaceAlpha',.2);
histogram(rev_rate_81Q,0:0.01:1,'Normalization','probability','FaceColor',colour_81Q,'EdgeColor',colour_81Q,'LineWidth',2,'FaceAlpha',.2);
xlabel('Reversal Rate'); ylabel('Frequency');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 1]); ylim([0 0.1]);

%Do t tests for velocity and run length distributions, edit below
[ant_vel_30Q_45Q,p_30Q_45Q_ant_vel]=ttest2(plus_vel_nonan_30Q,plus_vel_nonan_45Q);
[ant_vel_30Q_65Q,p_30Q_65Q_ant_vel]=ttest2(plus_vel_nonan_30Q,plus_vel_nonan_65Q);
[ant_vel_30Q_81Q,p_30Q_81Q_ant_vel]=ttest2(plus_vel_nonan_30Q,plus_vel_nonan_81Q);

[ret_vel_30Q_45Q,p_30Q_45Q_ret_vel]=ttest2(minus_vel_nonan_30Q,minus_vel_nonan_45Q);
[ret_vel_30Q_65Q,p_30Q_65Q_ret_vel]=ttest2(minus_vel_nonan_30Q,minus_vel_nonan_65Q);
[ret_vel_30Q_81Q,p_30Q_81Q_ret_vel]=ttest2(minus_vel_nonan_30Q,minus_vel_nonan_81Q);

[ant_lr_30Q_45Q,p_30Q_45Q_ant_lr]=ttest2(plus_lr_nonan_30Q,plus_lr_nonan_45Q);
[ant_lr_30Q_65Q,p_30Q_65Q_ant_lr]=ttest2(plus_lr_nonan_30Q,plus_lr_nonan_65Q);
[ant_lr_30Q_81Q,p_30Q_81Q_ant_lr]=ttest2(plus_lr_nonan_30Q,plus_lr_nonan_81Q);

[ret_lr_30Q_45Q,p_30Q_45Q_ret_lr]=ttest2(minus_lr_nonan_30Q,minus_lr_nonan_45Q);
[ret_lr_30Q_65Q,p_30Q_65Q_ret_lr]=ttest2(minus_lr_nonan_30Q,minus_lr_nonan_65Q);
[ret_lr_30Q_81Q,p_30Q_81Q_ret_lr]=ttest2(minus_lr_nonan_30Q,minus_lr_nonan_81Q);

titles={'30Q vs 30Q+IFg';'30Q vs 81Q';'30Q vs 81Q+IFg'};
% titles={'30Q vs 45Q';'30Q vs 65Q';'30Q vs 81Q'};
ant_vel_stats=[p_30Q_45Q_ant_vel;p_30Q_65Q_ant_vel;p_30Q_81Q_ant_vel];
ret_vel_stats=[p_30Q_45Q_ret_vel;p_30Q_65Q_ret_vel;p_30Q_81Q_ret_vel];
ant_lr_stats=[p_30Q_45Q_ant_lr;p_30Q_65Q_ant_lr;p_30Q_81Q_ant_lr];
ret_lr_stats=[p_30Q_45Q_ret_lr;p_30Q_65Q_ret_lr;p_30Q_81Q_ret_lr];

%ANOVA statistical testing
ant="Ant";
ret="Ret";

rl_ant_values=[plus_lr_nonan_30Q;plus_lr_nonan_45Q;plus_lr_nonan_65Q;plus_lr_nonan_81Q];
rl_ant_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(plus_lr_nonan_30Q) numel(plus_lr_nonan_45Q) numel(plus_lr_nonan_65Q) numel(plus_lr_nonan_81Q)])]';
rl_ant_data_tbl=table(rl_ant_values,rl_ant_groups);
[rl_ant_p,rl_ant_tbl_anova,rl_ant_stats]=anova1(rl_ant_values,rl_ant_groups,"off");
[rl_ant_results,a,b,gnames]=multcompare(rl_ant_stats);
rl_ant_tbl = array2table(rl_ant_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
rl_ant_tbl.("Group A") = gnames(rl_ant_tbl.("Group A"));
rl_ant_tbl.("Group B") = gnames(rl_ant_tbl.("Group B"));

rl_ret_values=[minus_lr_nonan_30Q;minus_lr_nonan_45Q;minus_lr_nonan_65Q;minus_lr_nonan_81Q];
rl_ret_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(minus_lr_nonan_30Q) numel(minus_lr_nonan_45Q) numel(minus_lr_nonan_65Q) numel(minus_lr_nonan_81Q)])]';
rl_ret_data_tbl=table(rl_ret_values,rl_ret_groups);
[rl_ret_p,rl_ret_tbl_anova,rl_ret_stats]=anova1(rl_ret_values,rl_ret_groups,"off");
[rl_ret_results,a,b,gnames]=multcompare(rl_ret_stats);
rl_ret_tbl = array2table(rl_ret_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
rl_ret_tbl.("Group A") = gnames(rl_ret_tbl.("Group A"));
rl_ret_tbl.("Group B") = gnames(rl_ret_tbl.("Group B"));

rl_values=[plus_lr_nonan_30Q;minus_lr_nonan_30Q;plus_lr_nonan_45Q;minus_lr_nonan_45Q;plus_lr_nonan_65Q;minus_lr_nonan_65Q;plus_lr_nonan_81Q;minus_lr_nonan_81Q];
cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(plus_lr_nonan_30Q)+numel(minus_lr_nonan_30Q) numel(plus_lr_nonan_45Q)+numel(minus_lr_nonan_45Q) numel(plus_lr_nonan_65Q)+numel(minus_lr_nonan_65Q) numel(plus_lr_nonan_81Q)+numel(minus_lr_nonan_81Q)])]';
directions=[repelem([ant, ret, ant, ret, ant, ret, ant, ret], [numel(plus_lr_nonan_30Q) numel(minus_lr_nonan_30Q) numel(plus_lr_nonan_45Q) numel(minus_lr_nonan_45Q) numel(plus_lr_nonan_65Q) numel(minus_lr_nonan_65Q) numel(plus_lr_nonan_81Q) numel(minus_lr_nonan_81Q)])]';
rl_data_tbl=table(rl_values,cell_types,directions);
[rl_p,rl_tbl_anova,rl_stats]=anovan(rl_values,{cell_types directions},"varnames",["cell type","direction"]);
[rl_results,a,b,gnames]=multcompare(rl_stats,"Dimension",[1 2],"Display","on");
rl_tbl = array2table(rl_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
rl_tbl.("Group A") = gnames(rl_tbl.("Group A"));
rl_tbl.("Group B") = gnames(rl_tbl.("Group B"));

vel_ant_values=[plus_vel_nonan_30Q;plus_vel_nonan_45Q;plus_vel_nonan_65Q;plus_vel_nonan_81Q];
vel_ant_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(plus_vel_nonan_30Q) numel(plus_vel_nonan_45Q) numel(plus_vel_nonan_65Q) numel(plus_vel_nonan_81Q)])]';
vel_ant_data_tbl=table(vel_ant_values,vel_ant_groups);
[vel_ant_p,vel_ant_tbl_anova,vel_ant_stats]=anova1(vel_ant_values,vel_ant_groups,"off");
[vel_ant_results,a,b,gnames]=multcompare(vel_ant_stats);
vel_ant_tbl = array2table(vel_ant_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
vel_ant_tbl.("Group A") = gnames(vel_ant_tbl.("Group A"));
vel_ant_tbl.("Group B") = gnames(vel_ant_tbl.("Group B"));

vel_ret_values=[minus_vel_nonan_30Q;minus_vel_nonan_45Q;minus_vel_nonan_65Q;minus_vel_nonan_81Q];
vel_ret_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(minus_vel_nonan_30Q) numel(minus_vel_nonan_45Q) numel(minus_vel_nonan_65Q) numel(minus_vel_nonan_81Q)])]';
vel_ret_data_tbl=table(vel_ret_values,vel_ret_groups);
[vel_ret_p,vel_ret_tbl_anova,vel_ret_stats]=anova1(vel_ret_values,vel_ret_groups,"off");
[vel_ret_results,a,b,gnames]=multcompare(vel_ret_stats);
vel_ret_tbl = array2table(vel_ret_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
vel_ret_tbl.("Group A") = gnames(vel_ret_tbl.("Group A"));
vel_ret_tbl.("Group B") = gnames(vel_ret_tbl.("Group B"));

prod_rl_vel_ant_values=[prod_rl_vel_30Q;prod_rl_vel_45Q;prod_rl_vel_65Q;prod_rl_vel_81Q];
prod_rl_vel_ant_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(prod_rl_vel_30Q) numel(prod_rl_vel_45Q) numel(prod_rl_vel_65Q) numel(prod_rl_vel_81Q)])]';
prod_rl_vel_ant_data_tbl=table(prod_rl_vel_ant_values,prod_rl_vel_ant_groups);
[prod_rl_vel_ant_p,prod_rl_vel_ant_tbl_anova,prod_rl_vel_ant_stats]=anova1(prod_rl_vel_ant_values,prod_rl_vel_ant_groups,"off");
[prod_rl_vel_ant_results,a,b,gnames]=multcompare(prod_rl_vel_ant_stats);
prod_rl_vel_ant_tbl = array2table(prod_rl_vel_ant_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
prod_rl_vel_ant_tbl.("Group A") = gnames(prod_rl_vel_ant_tbl.("Group A"));
prod_rl_vel_ant_tbl.("Group B") = gnames(prod_rl_vel_ant_tbl.("Group B"));

prod_rl_vel_ret_values=[prod_rl_vel_30Q_m;prod_rl_vel_45Q_m;prod_rl_vel_65Q_m;prod_rl_vel_81Q_m];
prod_rl_vel_ret_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(prod_rl_vel_30Q_m) numel(prod_rl_vel_45Q_m) numel(prod_rl_vel_65Q_m) numel(prod_rl_vel_81Q_m)])]';
prod_rl_vel_ret_data_tbl=table(prod_rl_vel_ret_values,prod_rl_vel_ant_groups);
[prod_rl_vel_ret_p,prod_rl_vel_ret_tbl_anova,prod_rl_vel_ret_stats]=anova1(prod_rl_vel_ret_values,prod_rl_vel_ret_groups,"off");
[prod_rl_vel_ret_results,a,b,gnames]=multcompare(prod_rl_vel_ret_stats);
prod_rl_vel_ret_tbl = array2table(prod_rl_vel_ret_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
prod_rl_vel_ret_tbl.("Group A") = gnames(prod_rl_vel_ret_tbl.("Group A"));
prod_rl_vel_ret_tbl.("Group B") = gnames(prod_rl_vel_ret_tbl.("Group B"));

vel_values=[plus_vel_nonan_30Q;minus_vel_nonan_30Q;plus_vel_nonan_45Q;minus_vel_nonan_45Q;plus_vel_nonan_65Q;minus_vel_nonan_65Q;plus_vel_nonan_81Q;minus_vel_nonan_81Q];
cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(plus_vel_nonan_30Q)+numel(minus_vel_nonan_30Q) numel(plus_vel_nonan_45Q)+numel(minus_vel_nonan_45Q) numel(plus_vel_nonan_65Q)+numel(minus_vel_nonan_65Q) numel(plus_vel_nonan_81Q)+numel(minus_vel_nonan_81Q)])]';
directions=[repelem([ant, ret, ant, ret, ant, ret, ant, ret], [numel(plus_vel_nonan_30Q) numel(minus_vel_nonan_30Q) numel(plus_vel_nonan_45Q) numel(minus_vel_nonan_45Q) numel(plus_vel_nonan_65Q) numel(minus_vel_nonan_65Q) numel(plus_vel_nonan_81Q) numel(minus_vel_nonan_81Q)])]';
vel_data_tbl=table(vel_values,cell_types,directions);
[vel_p,vel_tbl_anova,vel_stats]=anovan(vel_values,{cell_types directions},"varnames",["cell type","direction"]);
[vel_results,a,b,gnames]=multcompare(vel_stats,"Dimension",[1 2],"Display","on");
vel_tbl = array2table(vel_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
vel_tbl.("Group A") = gnames(vel_tbl.("Group A"));
vel_tbl.("Group B") = gnames(vel_tbl.("Group B"));

fileprefix='20240425_isoHD_KB_mito_per_traj_';
% 
tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
FolderName_1 = tempdir_1;   % Your destination folder
% %Saving stats to a table
writetable(rl_ant_tbl,fullfile(FolderName_1, [fileprefix,'rl_ant.csv']),'WriteRowNames',true);
writetable(rl_ant_data_tbl,fullfile(FolderName_1, [fileprefix,'rl_ant_data.csv']),'WriteRowNames',true);
writetable(rl_ret_tbl,fullfile(FolderName_1, [fileprefix,'rl_ret.csv']),'WriteRowNames',true);
writetable(rl_ret_data_tbl,fullfile(FolderName_1, [fileprefix,'rl_ret_data.csv']),'WriteRowNames',true);
writetable(rl_tbl,fullfile(FolderName_1, [fileprefix,'rl.csv']),'WriteRowNames',true);
writetable(rl_data_tbl,fullfile(FolderName_1,[fileprefix,'rl_data.csv']),'WriteRowNames',true);
writetable(vel_ant_tbl,fullfile(FolderName_1, [fileprefix,'vel_ant.csv']),'WriteRowNames',true);
writetable(vel_ant_data_tbl, fullfile(FolderName_1, [fileprefix,'vel_ant_data.csv']),'WriteRowNames',true);
writetable(vel_ret_tbl,fullfile(FolderName_1, [fileprefix,'vel_ret.csv']),'WriteRowNames',true);
writetable(vel_ret_data_tbl, fullfile(FolderName_1, [fileprefix,'vel_ret_data.csv']),'WriteRowNames',true);
writetable(prod_rl_vel_ant_tbl,fullfile(FolderName_1, [fileprefix,'prod_rl_vel_ant.csv']),'WriteRowNames',true);
writetable(prod_rl_vel_ant_data_tbl, fullfile(FolderName_1, [fileprefix,'prod_rl_vel_ant_data.csv']),'WriteRowNames',true);
writetable(prod_rl_vel_ret_tbl,fullfile(FolderName_1, [fileprefix,'prod_rl_vel_ret.csv']),'WriteRowNames',true);
writetable(prod_rl_vel_ret_data_tbl, fullfile(FolderName_1, [fileprefix,'prod_rl_vel_ret_data.csv']),'WriteRowNames',true);
writetable(vel_tbl,fullfile(FolderName_1, [fileprefix,'vel.csv']),'WriteRowNames',true);
writetable(vel_data_tbl, fullfile(FolderName_1, [fileprefix,'vel_data.csv']),'WriteRowNames',true);

% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end