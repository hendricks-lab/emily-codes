
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');

close all;

% polyQ length
% cmap = colormap(hot(100));
% colour_30Q=[0.5 0.5 0.5];
% colour_45Q=cmap(20,:);%[0.75 0 0];
% colour_65Q=cmap(25,:);%[0.5 0 0];
% colour_81Q=cmap(30,:);%[0.25 0 0];
% htt_30Q="30Q";
% htt_45Q="45Q";
% htt_65Q="65Q";
% htt_81Q="81Q";
% x_labels={'30Q','45Q', '65Q','81Q'};

%stress condition
cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_65Q=cmap(30,:);%[0.5 0 0];
cmap2= colormap(spring(100));
colour_45Q=[1.0000 0.4444 0.5556];%cmap2(50,:);%[0.75 0 0];
colour_81Q=cmap2(30,:);%[0.25 0 0];
htt_30Q="30Q";
htt_45Q="30Q+IFg";
htt_65Q="81Q";
htt_81Q="81Q+IFg";
x_labels={'30Q','30Q+IFNγ', '81Q','81Q+IFNγ'};

% pth='/Volumes/Emily_htt_2/Neuron/bdnf_flux/';
pth='/Volumes/Emily_htt_2/Neuron/lyso_flux/';
% pth='/Volumes/Emily_htt_2/Neuron/mito_flux/';

fileprefix='20240119_isoHD_lyso_IFg_';

tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/';  % Your destination folder for figures
tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder for stats

load(fullfile(pth,'KB_frac_ant_ex_IFg_5um.mat'));
load(fullfile(pth,'KB_frac_ret_ex_IFg_5um.mat'));
load(fullfile(pth,'KB_total_enters_IFg_5um.mat'));
load(fullfile(pth,'KB_frac_stay_in_IFg_5um.mat'));


[ant_30Q_45Q,p_30Q_45Q_ant]=ttest2(frac_ant_exits_per_cell_all_conds{1},frac_ant_exits_per_cell_all_conds{2});
[ant_30Q_65Q,p_30Q_65Q_ant]=ttest2(frac_ant_exits_per_cell_all_conds{1},frac_ant_exits_per_cell_all_conds{3});
[ant_30Q_81Q,p_30Q_81Q_ant]=ttest2(frac_ant_exits_per_cell_all_conds{1},frac_ant_exits_per_cell_all_conds{4});

[ret_30Q_45Q,p_30Q_45Q_ret]=ttest2(frac_ret_exits_per_cell_all_conds{1},frac_ret_exits_per_cell_all_conds{2});
[ret_30Q_65Q,p_30Q_65Q_ret]=ttest2(frac_ret_exits_per_cell_all_conds{1},frac_ret_exits_per_cell_all_conds{3});
[ret_30Q_81Q,p_30Q_81Q_ret]=ttest2(frac_ret_exits_per_cell_all_conds{1},frac_ret_exits_per_cell_all_conds{4});

[nonet_30Q_45Q,p_30Q_45Q_nonet]=ttest2(frac_stays_in_per_cell_all_conds{1},frac_stays_in_per_cell_all_conds{2});
[nonet_30Q_65Q,p_30Q_65Q_nonet]=ttest2(frac_stays_in_per_cell_all_conds{1},frac_stays_in_per_cell_all_conds{3});
[nonet_30Q_81Q,p_30Q_81Q_nonet]=ttest2(frac_stays_in_per_cell_all_conds{1},frac_stays_in_per_cell_all_conds{4});

[tot_ent_30Q_45Q,p_30Q_45Q_tot_ent]=ttest2(tot_enters_per_cell_all_conds{1},tot_enters_per_cell_all_conds{2});
[tot_ent_30Q_65Q,p_30Q_65Q_tot_ent]=ttest2(tot_enters_per_cell_all_conds{1},tot_enters_per_cell_all_conds{3});
[tot_ent_30Q_81Q,p_30Q_81Q_tot_ent]=ttest2(tot_enters_per_cell_all_conds{1},tot_enters_per_cell_all_conds{4});

% titles={'30Q vs 30Q+IFg';'30Q vs 81Q';'30Q vs 81Q+IFg'};
% % titles={'30Q vs 45Q';'30Q vs 65Q';'30Q vs 81Q'};
% ant_stats=[p_30Q_45Q_ant;p_30Q_65Q_ant;p_30Q_81Q_ant];
% ret_stats=[p_30Q_45Q_ret;p_30Q_65Q_ret;p_30Q_81Q_ret];
% nonet_stats=[p_30Q_45Q_nonet;p_30Q_65Q_nonet;p_30Q_81Q_nonet];
% tot_ent_stats=[p_30Q_45Q_tot_ent;p_30Q_65Q_tot_ent;p_30Q_81Q_tot_ent];
% 
% 
% FolderName_1 = tempdir_1;   % Your destination folder
% %Saving stats to a table
% stats_table=table(titles, ant_stats,ret_stats, nonet_stats, tot_ent_stats);
% writetable(stats_table,fullfile(FolderName_1, [fileprefix,'.csv']),'WriteRowNames',true);

ant_vs_ret_raw{1}=frac_ant_exits_per_cell_all_conds{1}./frac_ret_exits_per_cell_all_conds{1};
ant_vs_ret_raw{2}=frac_ant_exits_per_cell_all_conds{2}./frac_ret_exits_per_cell_all_conds{2};
ant_vs_ret_raw{3}=frac_ant_exits_per_cell_all_conds{3}./frac_ret_exits_per_cell_all_conds{3};
ant_vs_ret_raw{4}=frac_ant_exits_per_cell_all_conds{4}./frac_ret_exits_per_cell_all_conds{4};

ant_vs_ret{1}=ant_vs_ret_raw{1}(isfinite(ant_vs_ret_raw{1}));
ant_vs_ret{2}=ant_vs_ret_raw{2}(isfinite(ant_vs_ret_raw{2}));
ant_vs_ret{3}=ant_vs_ret_raw{3}(isfinite(ant_vs_ret_raw{3}));
ant_vs_ret{4}=ant_vs_ret_raw{4}(isfinite(ant_vs_ret_raw{4}));

%Calculating confidence intervals
SEM_tot_ent{1}=std(tot_enters_per_cell_all_conds{1})/sqrt(length(tot_enters_per_cell_all_conds{1}));
SEM_tot_ent{2}=std(tot_enters_per_cell_all_conds{2})/sqrt(length(tot_enters_per_cell_all_conds{2}));
SEM_tot_ent{3}=std(tot_enters_per_cell_all_conds{3})/sqrt(length(tot_enters_per_cell_all_conds{3}));
SEM_tot_ent{4}=std(tot_enters_per_cell_all_conds{4})/sqrt(length(tot_enters_per_cell_all_conds{4}));

SEM_tot_ent_pm{1}=std(tot_enters_per_cell_all_conds{1}./25)/sqrt(length(tot_enters_per_cell_all_conds{1}));
SEM_tot_ent_pm{2}=std(tot_enters_per_cell_all_conds{2}./25)/sqrt(length(tot_enters_per_cell_all_conds{2}));
SEM_tot_ent_pm{3}=std(tot_enters_per_cell_all_conds{3}./25)/sqrt(length(tot_enters_per_cell_all_conds{3}));
SEM_tot_ent_pm{4}=std(tot_enters_per_cell_all_conds{4}./25)/sqrt(length(tot_enters_per_cell_all_conds{4}));

SEM_ant_ret{1}=std(ant_vs_ret{1})/sqrt(length(ant_vs_ret{1}));
SEM_ant_ret{2}=std(ant_vs_ret{2})/sqrt(length(ant_vs_ret{2}));
SEM_ant_ret{3}=std(ant_vs_ret{3})/sqrt(length(ant_vs_ret{3}));
SEM_ant_ret{4}=std(ant_vs_ret{4})/sqrt(length(ant_vs_ret{4}));

ts_tot_ent{1} = tinv([0.025  0.975],length(tot_enters_per_cell_all_conds{1})-1);      % T-Score
ts_tot_ent{2} = tinv([0.025  0.975],length(tot_enters_per_cell_all_conds{2})-1); 
ts_tot_ent{3} = tinv([0.025  0.975],length(tot_enters_per_cell_all_conds{3})-1); 
ts_tot_ent{4} = tinv([0.025  0.975],length(tot_enters_per_cell_all_conds{4})-1); 

ts_ant_ret{1} = tinv([0.025  0.975],length(ant_vs_ret{1})-1);      % T-Score
ts_ant_ret{2} = tinv([0.025  0.975],length(ant_vs_ret{2})-1); 
ts_ant_ret{3} = tinv([0.025  0.975],length(ant_vs_ret{3})-1); 
ts_ant_ret{4} = tinv([0.025  0.975],length(ant_vs_ret{4})-1); 

CI_tot_ent{1} = mean(tot_enters_per_cell_all_conds{1}) + ts_tot_ent{1}*SEM_tot_ent{1};   
CI_tot_ent{2} = mean(tot_enters_per_cell_all_conds{2}) + ts_tot_ent{2}*SEM_tot_ent{2};
CI_tot_ent{3} = mean(tot_enters_per_cell_all_conds{3}) + ts_tot_ent{3}*SEM_tot_ent{3}; 
CI_tot_ent{4} = mean(tot_enters_per_cell_all_conds{4}) + ts_tot_ent{4}*SEM_tot_ent{4}; 

CI_ant_ret{1} = mean(ant_vs_ret{1}) + ts_ant_ret{1}*SEM_ant_ret{1};   
CI_ant_ret{2} = mean(ant_vs_ret{2}) + ts_ant_ret{2}*SEM_ant_ret{2};
CI_ant_ret{3} = mean(ant_vs_ret{3}) + ts_ant_ret{3}*SEM_ant_ret{3}; 
CI_ant_ret{4} = mean(ant_vs_ret{4}) + ts_ant_ret{4}*SEM_ant_ret{4}; 

%Setting a threshold for maximum total cargoes
tot_enters_thres=100;
tot_ent_b_thres{1}=tot_enters_per_cell_all_conds{1}(find(tot_enters_per_cell_all_conds{1}<=tot_enters_thres));
tot_ent_b_thres{2}=tot_enters_per_cell_all_conds{2}(find(tot_enters_per_cell_all_conds{2}<=tot_enters_thres));
tot_ent_b_thres{3}=tot_enters_per_cell_all_conds{3}(find(tot_enters_per_cell_all_conds{3}<=tot_enters_thres));
tot_ent_b_thres{4}=tot_enters_per_cell_all_conds{4}(find(tot_enters_per_cell_all_conds{4}<=tot_enters_thres));

tot_ent_a_thres{1}=tot_enters_per_cell_all_conds{1}(find(tot_enters_per_cell_all_conds{1}>tot_enters_thres));
tot_ent_a_thres{2}=tot_enters_per_cell_all_conds{2}(find(tot_enters_per_cell_all_conds{2}>tot_enters_thres));
tot_ent_a_thres{3}=tot_enters_per_cell_all_conds{3}(find(tot_enters_per_cell_all_conds{3}>tot_enters_thres));
tot_ent_a_thres{4}=tot_enters_per_cell_all_conds{4}(find(tot_enters_per_cell_all_conds{4}>tot_enters_thres));

tot_enters_thres_pm=4;
tot_enters_per_cell_all_conds_pm{1}=tot_enters_per_cell_all_conds{1}./25;
tot_enters_per_cell_all_conds_pm{2}=tot_enters_per_cell_all_conds{2}./25;
tot_enters_per_cell_all_conds_pm{3}=tot_enters_per_cell_all_conds{3}./25;
tot_enters_per_cell_all_conds_pm{4}=tot_enters_per_cell_all_conds{4}./25;

tot_ent_b_thres_pm{1}=tot_enters_per_cell_all_conds_pm{1}(find(tot_enters_per_cell_all_conds_pm{1}<tot_enters_thres_pm));
tot_ent_b_thres_pm{2}=tot_enters_per_cell_all_conds_pm{2}(find(tot_enters_per_cell_all_conds_pm{2}<tot_enters_thres_pm));
tot_ent_b_thres_pm{3}=tot_enters_per_cell_all_conds_pm{3}(find(tot_enters_per_cell_all_conds_pm{3}<tot_enters_thres_pm));
tot_ent_b_thres_pm{4}=tot_enters_per_cell_all_conds_pm{4}(find(tot_enters_per_cell_all_conds_pm{4}<tot_enters_thres_pm));

tot_ent_a_thres_pm{1}=tot_enters_per_cell_all_conds_pm{1}(find(tot_enters_per_cell_all_conds_pm{1}>=tot_enters_thres_pm));
tot_ent_a_thres_pm{2}=tot_enters_per_cell_all_conds_pm{2}(find(tot_enters_per_cell_all_conds_pm{2}>=tot_enters_thres_pm));
tot_ent_a_thres_pm{3}=tot_enters_per_cell_all_conds_pm{3}(find(tot_enters_per_cell_all_conds_pm{3}>=tot_enters_thres_pm));
tot_ent_a_thres_pm{4}=tot_enters_per_cell_all_conds_pm{4}(find(tot_enters_per_cell_all_conds_pm{4}>=tot_enters_thres_pm));

%Setting a threshold for maximum ant/ret ratio
ant_ret_thres=5;
ant_ret_b_thres{1}=ant_vs_ret{1}(find(ant_vs_ret{1}<ant_ret_thres));
ant_ret_b_thres{2}=ant_vs_ret{2}(find(ant_vs_ret{2}<ant_ret_thres));
ant_ret_b_thres{3}=ant_vs_ret{3}(find(ant_vs_ret{3}<ant_ret_thres));
ant_ret_b_thres{4}=ant_vs_ret{4}(find(ant_vs_ret{4}<ant_ret_thres));

ant_ret_a_thres{1}=ant_vs_ret{1}(find(ant_vs_ret{1}>=ant_ret_thres));
ant_ret_a_thres{2}=ant_vs_ret{2}(find(ant_vs_ret{2}>=ant_ret_thres));
ant_ret_a_thres{3}=ant_vs_ret{3}(find(ant_vs_ret{3}>=ant_ret_thres));
ant_ret_a_thres{4}=ant_vs_ret{4}(find(ant_vs_ret{4}>=ant_ret_thres));

%Making variables easier for bar plots
avg_total_entered(1)=mean(tot_enters_per_cell_all_conds{1},"omitnan");
avg_total_entered(2)=mean(tot_enters_per_cell_all_conds{2},"omitnan");
avg_total_entered(3)=mean(tot_enters_per_cell_all_conds{3},"omitnan");
avg_total_entered(4)=mean(tot_enters_per_cell_all_conds{4},"omitnan");
avg_exits_ant(1)=mean(frac_ant_exits_per_cell_all_conds{1},"omitnan");
avg_exits_ant(2)=mean(frac_ant_exits_per_cell_all_conds{2},"omitnan");
avg_exits_ant(3)=mean(frac_ant_exits_per_cell_all_conds{3},"omitnan");
avg_exits_ant(4)=mean(frac_ant_exits_per_cell_all_conds{4},"omitnan");
avg_exits_ret(1)=mean(frac_ret_exits_per_cell_all_conds{1},"omitnan");
avg_exits_ret(2)=mean(frac_ret_exits_per_cell_all_conds{2},"omitnan");
avg_exits_ret(3)=mean(frac_ret_exits_per_cell_all_conds{3},"omitnan");
avg_exits_ret(4)=mean(frac_ret_exits_per_cell_all_conds{4},"omitnan");
avg_stays_in(1)=mean(frac_stays_in_per_cell_all_conds{1},"omitnan");
avg_stays_in(2)=mean(frac_stays_in_per_cell_all_conds{2},"omitnan");
avg_stays_in(3)=mean(frac_stays_in_per_cell_all_conds{3},"omitnan");
avg_stays_in(4)=mean(frac_stays_in_per_cell_all_conds{4},"omitnan");

sem_total_entered(1)=std(tot_enters_per_cell_all_conds{1},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{1})));
sem_total_entered(2)=std(tot_enters_per_cell_all_conds{2},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{2})));
sem_total_entered(3)=std(tot_enters_per_cell_all_conds{3},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{3})));
sem_total_entered(4)=std(tot_enters_per_cell_all_conds{4},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{4})));

sem_total_entered_pm(1)=std(tot_enters_per_cell_all_conds_pm{1},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{1})));
sem_total_entered_pm(2)=std(tot_enters_per_cell_all_conds_pm{2},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{2})));
sem_total_entered_pm(3)=std(tot_enters_per_cell_all_conds_pm{3},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{3})));
sem_total_entered_pmm(4)=std(tot_enters_per_cell_all_conds_pm{4},"omitnan")/(sqrt(length(tot_enters_per_cell_all_conds{4})));

sem_exits_ant(1)=std(frac_ant_exits_per_cell_all_conds{1},"omitnan")/(sqrt(length(frac_ant_exits_per_cell_all_conds{1})));
sem_exits_ant(2)=std(frac_ant_exits_per_cell_all_conds{2},"omitnan")/(sqrt(length(frac_ant_exits_per_cell_all_conds{2})));
sem_exits_ant(3)=std(frac_ant_exits_per_cell_all_conds{3},"omitnan")/(sqrt(length(frac_ant_exits_per_cell_all_conds{3})));
sem_exits_ant(4)=std(frac_ant_exits_per_cell_all_conds{4},"omitnan")/(sqrt(length(frac_ant_exits_per_cell_all_conds{4})));
sem_exits_ret(1)=std(frac_ret_exits_per_cell_all_conds{1},"omitnan")/(sqrt(length(frac_ret_exits_per_cell_all_conds{1})));
sem_exits_ret(2)=std(frac_ret_exits_per_cell_all_conds{2},"omitnan")/(sqrt(length(frac_ret_exits_per_cell_all_conds{2})));
sem_exits_ret(3)=std(frac_ret_exits_per_cell_all_conds{3},"omitnan")/(sqrt(length(frac_ret_exits_per_cell_all_conds{3})));
sem_exits_ret(4)=std(frac_ret_exits_per_cell_all_conds{4},"omitnan")/(sqrt(length(frac_ret_exits_per_cell_all_conds{4})));
sem_stays_in(1)=std(frac_stays_in_per_cell_all_conds{1},"omitnan")/(sqrt(length(frac_stays_in_per_cell_all_conds{1})));
sem_stays_in(2)=std(frac_stays_in_per_cell_all_conds{2},"omitnan")/(sqrt(length(frac_stays_in_per_cell_all_conds{2})));
sem_stays_in(3)=std(frac_stays_in_per_cell_all_conds{3},"omitnan")/(sqrt(length(frac_stays_in_per_cell_all_conds{3})));
sem_stays_in(4)=std(frac_stays_in_per_cell_all_conds{4},"omitnan")/(sqrt(length(frac_stays_in_per_cell_all_conds{4})));

figure('Name','ant_ret_stat_exit_bar','NumberTitle','off'), hold on, %4 conditions
exit_fractions=[avg_exits_ant(1) avg_exits_ant(2) avg_exits_ant(3) avg_exits_ant(4); avg_exits_ret(1) avg_exits_ret(2) avg_exits_ret(3) avg_exits_ret(4); avg_stays_in(1) avg_stays_in(2) avg_stays_in(3) avg_stays_in(4)];
exit_err_bars=[sem_exits_ant(1) sem_exits_ant(2) sem_exits_ant(3) sem_exits_ant(4); sem_exits_ret(1) sem_exits_ret(2) sem_exits_ret(3) sem_exits_ret(4); sem_stays_in(1) sem_stays_in(2) sem_stays_in(3) sem_stays_in(4)];
b1=bar(exit_fractions,'grouped');
b1(1).FaceColor=colour_30Q;
b1(2).FaceColor=colour_45Q;
b1(3).FaceColor=colour_65Q;
b1(4).FaceColor=colour_81Q;
b1(1).FaceAlpha=0.5;
b1(2).FaceAlpha=0.5;
b1(3).FaceAlpha=0.5;
b1(4).FaceAlpha=0.5;
hold on
scatter(0.7273-0.02+0.05*rand(numel(frac_ant_exits_per_cell_all_conds{1}),1),frac_ant_exits_per_cell_all_conds{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(0.9091-0.02+0.05*rand(numel(frac_ant_exits_per_cell_all_conds{2}),1),frac_ant_exits_per_cell_all_conds{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.0909-0.02+0.05*rand(numel(frac_ant_exits_per_cell_all_conds{3}),1),frac_ant_exits_per_cell_all_conds{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.2727-0.02+0.05*rand(numel(frac_ant_exits_per_cell_all_conds{4}),1),frac_ant_exits_per_cell_all_conds{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.7273-0.02+0.05*rand(numel(frac_ret_exits_per_cell_all_conds{1}),1),frac_ret_exits_per_cell_all_conds{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.9091-0.02+0.05*rand(numel(frac_ret_exits_per_cell_all_conds{2}),1),frac_ret_exits_per_cell_all_conds{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.0909-0.02+0.05*rand(numel(frac_ret_exits_per_cell_all_conds{3}),1),frac_ret_exits_per_cell_all_conds{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.2727-0.02+0.05*rand(numel(frac_ret_exits_per_cell_all_conds{4}),1),frac_ret_exits_per_cell_all_conds{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.7273-0.02+0.05*rand(numel(frac_stays_in_per_cell_all_conds{1}),1),frac_stays_in_per_cell_all_conds{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.9091-0.02+0.05*rand(numel(frac_stays_in_per_cell_all_conds{2}),1),frac_stays_in_per_cell_all_conds{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3.0909-0.02+0.05*rand(numel(frac_stays_in_per_cell_all_conds{3}),1),frac_stays_in_per_cell_all_conds{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3.2727-0.02+0.05*rand(numel(frac_stays_in_per_cell_all_conds{4}),1),frac_stays_in_per_cell_all_conds{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
hold on
[ngroups,nbars]= size(exit_fractions);
groupwidth= min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
   % Calculate center of each bar
   x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
   errorbar(x, exit_fractions (:,i), exit_err_bars (:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
xticks([1 2 3]);
xticklabels({'Ant flux', 'Ret flux','No net flux'});
ylabel('Fraction of Cargoes');
publication_fig(0,0,1);
set(gca, 'FontSize', 40);
box("off");

figure('Name','total_enters','NumberTitle','off'), hold on,
plot(1,mean(tot_enters_per_cell_all_conds{1}),'k_','MarkerSize',30);
plot(2,mean(tot_enters_per_cell_all_conds{2}),'k_','MarkerSize',30);
plot(3,mean(tot_enters_per_cell_all_conds{3}),'k_','MarkerSize',30);
plot(4,mean(tot_enters_per_cell_all_conds{4}),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(tot_ent_b_thres{1}),1),tot_ent_b_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(tot_ent_b_thres{2}),1),tot_ent_b_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(tot_ent_b_thres{3}),1),tot_ent_b_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(tot_ent_b_thres{4}),1),tot_ent_b_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(tot_ent_a_thres{1}),1),(tot_enters_thres+25)*(ones(numel(tot_ent_a_thres{1}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(tot_ent_a_thres{2}),1),(tot_enters_thres+25)*(ones(numel(tot_ent_a_thres{2}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(tot_ent_a_thres{3}),1),(tot_enters_thres+25)*(ones(numel(tot_ent_a_thres{3}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(tot_ent_a_thres{4}),1),(tot_enters_thres+25)*(ones(numel(tot_ent_a_thres{4}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(tot_enters_per_cell_all_conds{1}),SEM_tot_ent{1},'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(tot_enters_per_cell_all_conds{2}),SEM_tot_ent{2},'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(tot_enters_per_cell_all_conds{3}),SEM_tot_ent{3},'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(tot_enters_per_cell_all_conds{4}),SEM_tot_ent{4},'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 5]);
xticks([1 2 3 4]);
xticklabels(x_labels);
ylabel('Total Cargoes');
publication_fig(0,0,1);
box("off");
pbaspect([1 1 1]);
set(gca, 'FontSize',40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','total_enters_pm','NumberTitle','off'), hold on,
plot(1,mean(tot_enters_per_cell_all_conds_pm{1}),'k_','MarkerSize',30);
plot(2,mean(tot_enters_per_cell_all_conds_pm{2}),'k_','MarkerSize',30);
plot(3,mean(tot_enters_per_cell_all_conds_pm{3}),'k_','MarkerSize',30);
plot(4,mean(tot_enters_per_cell_all_conds_pm{4}),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(tot_ent_b_thres_pm{1}),1),tot_ent_b_thres_pm{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(tot_ent_b_thres_pm{2}),1),tot_ent_b_thres_pm{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(tot_ent_b_thres_pm{3}),1),tot_ent_b_thres_pm{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(tot_ent_b_thres_pm{4}),1),tot_ent_b_thres_pm{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(tot_ent_a_thres_pm{1}),1),(tot_enters_thres_pm+0.2)*(ones(numel(tot_ent_a_thres_pm{1}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(tot_ent_a_thres_pm{2}),1),(tot_enters_thres_pm+0.2)*(ones(numel(tot_ent_a_thres_pm{2}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(tot_ent_a_thres_pm{3}),1),(tot_enters_thres_pm+0.2)*(ones(numel(tot_ent_a_thres_pm{3}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(tot_ent_a_thres_pm{4}),1),(tot_enters_thres_pm+0.2)*(ones(numel(tot_ent_a_thres_pm{4}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(tot_enters_per_cell_all_conds_pm{1}),SEM_tot_ent_pm{1},'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(tot_enters_per_cell_all_conds_pm{2}),SEM_tot_ent_pm{2},'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(tot_enters_per_cell_all_conds_pm{3}),SEM_tot_ent_pm{3},'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(tot_enters_per_cell_all_conds_pm{4}),SEM_tot_ent_pm{4},'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 4.5]);
xticks([1 2 3 4]);
xticklabels(x_labels);
ylabel('Cargoes/\mum');
publication_fig(0,0,1);
box("off");
pbaspect([1 1 1]);
set(gca, 'FontSize',40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');


figure('Name','ant_ret_thres','NumberTitle','off'), hold on,
plot(1,mean(ant_vs_ret{1}),'k_','MarkerSize',30);
plot(2,mean(ant_vs_ret{2}),'k_','MarkerSize',30);
plot(3,mean(ant_vs_ret{3}),'k_','MarkerSize',30);
plot(4,mean(ant_vs_ret{4}),'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(ant_ret_b_thres{1}),1),ant_ret_b_thres{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(ant_ret_b_thres{2}),1),ant_ret_b_thres{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(ant_ret_b_thres{3}),1),ant_ret_b_thres{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(ant_ret_b_thres{4}),1),ant_ret_b_thres{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(ant_ret_a_thres{1}),1),(ant_ret_thres+0.2)*(ones(numel(ant_ret_a_thres{1}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(ant_ret_a_thres{2}),1),(ant_ret_thres+0.2)*(ones(numel(ant_ret_a_thres{2}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(ant_ret_a_thres{3}),1),(ant_ret_thres+0.2)*(ones(numel(ant_ret_a_thres{3}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(ant_ret_a_thres{4}),1),(ant_ret_thres+0.2)*(ones(numel(ant_ret_a_thres{4}),1)),50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(ant_vs_ret{1}),SEM_ant_ret{1},'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(ant_vs_ret{2}),SEM_ant_ret{2},'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(ant_vs_ret{3}),SEM_ant_ret{3},'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(ant_vs_ret{4}),SEM_ant_ret{4},'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 5]);
xticks([1 2 3 4]);
xticklabels(x_labels);
ylabel('Ant/Ret');
publication_fig(0,0,1);
box("off");
pbaspect([1 1 1]);
set(gca, 'FontSize',40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','ant_vs_ret','NumberTitle','off'), hold on,
plot(1,mean(ant_vs_ret{1}),'k_','MarkerSize',30);
plot(2,mean(ant_vs_ret{2}),'k_','MarkerSize',30);
plot(3,mean(ant_vs_ret{3}),'k_','MarkerSize',30);
plot(4,mean(ant_vs_ret{4}),'k+','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(ant_vs_ret{1}),1),ant_vs_ret{1},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(ant_vs_ret{2}),1),ant_vs_ret{2},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(ant_vs_ret{3}),1),ant_vs_ret{3},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(ant_vs_ret{4}),1),ant_vs_ret{4},50,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
errorbar(1,mean(ant_vs_ret{1}),SEM_ant_ret{1},'Color',colour_30Q,'LineWidth',4,'CapSize',0);
errorbar(2,mean(ant_vs_ret{2}),SEM_ant_ret{2},'Color',colour_45Q,'LineWidth',4,'CapSize',0);
errorbar(3,mean(ant_vs_ret{3}),SEM_ant_ret{3},'Color',colour_65Q,'LineWidth',4,'CapSize',0);
errorbar(4,mean(ant_vs_ret{4}),SEM_ant_ret{4},'Color',colour_81Q,'LineWidth',4,'CapSize',0);
xlim([0 5]);
xticks([1 2 3 4]);
xticklabels(x_labels);
ylabel('Ant/Ret');
publication_fig(0,0,1);
box("off");
set(gca, 'FontSize', 40);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

%ANOVAs for all 5um flux analysis
ant="Ant";
ret="Ret";
no_net="No net";

frac_values=[frac_ant_exits_per_cell_all_conds{1}';frac_ret_exits_per_cell_all_conds{1}';frac_stays_in_per_cell_all_conds{1}'; frac_ant_exits_per_cell_all_conds{2}';frac_ret_exits_per_cell_all_conds{2}';frac_stays_in_per_cell_all_conds{2}';frac_ant_exits_per_cell_all_conds{3}';frac_ret_exits_per_cell_all_conds{3}';frac_stays_in_per_cell_all_conds{3}';frac_ant_exits_per_cell_all_conds{4}';frac_ret_exits_per_cell_all_conds{4}';frac_stays_in_per_cell_all_conds{4}'];
cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(frac_ant_exits_per_cell_all_conds{1})+numel(frac_ret_exits_per_cell_all_conds{1})+numel(frac_stays_in_per_cell_all_conds{1}) numel(frac_ant_exits_per_cell_all_conds{2})+numel(frac_ret_exits_per_cell_all_conds{2})+numel(frac_stays_in_per_cell_all_conds{2}) numel(frac_ant_exits_per_cell_all_conds{3})+numel(frac_ret_exits_per_cell_all_conds{3})+numel(frac_stays_in_per_cell_all_conds{3}) numel(frac_ant_exits_per_cell_all_conds{4})+numel(frac_ret_exits_per_cell_all_conds{4})+numel(frac_stays_in_per_cell_all_conds{4})])]';
flux_types=[repelem([ant, ret,no_net, ant, ret, no_net,ant, ret, no_net,ant, ret, no_net], [numel(frac_ant_exits_per_cell_all_conds{1}) numel(frac_ret_exits_per_cell_all_conds{1}) numel(frac_stays_in_per_cell_all_conds{1}) numel(frac_ant_exits_per_cell_all_conds{2}) numel(frac_ret_exits_per_cell_all_conds{2}) numel(frac_stays_in_per_cell_all_conds{2}) numel(frac_ant_exits_per_cell_all_conds{3}) numel(frac_ret_exits_per_cell_all_conds{3}) numel(frac_stays_in_per_cell_all_conds{3}) numel(frac_ant_exits_per_cell_all_conds{4}) numel(frac_ret_exits_per_cell_all_conds{4}) numel(frac_stays_in_per_cell_all_conds{4})])]';
frac_data_tbl=table(frac_values,cell_types,flux_types);
[frac_p,frac_tbl_anova,frac_stats]=anovan(frac_values,{cell_types flux_types},"varnames",["cell type","direction"]);
[frac_results,a,b,gnames]=multcompare(frac_stats,"Dimension",[1 2],"Display","off");
frac_tbl = array2table(frac_results,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
frac_tbl.("Group A")=gnames(frac_tbl.("Group A"));
frac_tbl.("Group B")=gnames(frac_tbl.("Group B"));

ant_frac_values=[frac_ant_exits_per_cell_all_conds{1}'; frac_ant_exits_per_cell_all_conds{2}';frac_ant_exits_per_cell_all_conds{3}';frac_ant_exits_per_cell_all_conds{4}'];
ant_frac_cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(frac_ant_exits_per_cell_all_conds{1}) numel(frac_ant_exits_per_cell_all_conds{2}) numel(frac_ant_exits_per_cell_all_conds{3}) numel(frac_ant_exits_per_cell_all_conds{4})])]';
ant_frac_data_tbl=table(ant_frac_values,ant_frac_cell_types);
[ant_frac_p,ant_frac_tbl_anova,ant_frac_stats]=anova1(ant_frac_values,ant_frac_cell_types,"on");
[ant_frac_results,a,b,gnames]=multcompare(ant_frac_stats);
ant_frac_tbl = array2table(ant_frac_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
ant_frac_tbl.("Group A")=gnames(ant_frac_tbl.("Group A"));
ant_frac_tbl.("Group B")=gnames(ant_frac_tbl.("Group B"));

ret_frac_values=[frac_ret_exits_per_cell_all_conds{1}'; frac_ret_exits_per_cell_all_conds{2}';frac_ret_exits_per_cell_all_conds{3}';frac_ret_exits_per_cell_all_conds{4}'];
ret_frac_cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(frac_ret_exits_per_cell_all_conds{1}) numel(frac_ret_exits_per_cell_all_conds{2}) numel(frac_ret_exits_per_cell_all_conds{3}) numel(frac_ret_exits_per_cell_all_conds{4})])]';
ret_frac_data_tbl=table(ret_frac_values,ret_frac_cell_types);
[ret_frac_p,ret_frac_tbl_anova,ret_frac_stats]=anova1(ret_frac_values,ret_frac_cell_types,"on");
[ret_frac_results,a,b,gnames]=multcompare(ret_frac_stats);
ret_frac_tbl = array2table(ret_frac_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
ret_frac_tbl.("Group A")=gnames(ret_frac_tbl.("Group A"));
ret_frac_tbl.("Group B")=gnames(ret_frac_tbl.("Group B"));

no_net_frac_values=[frac_stays_in_per_cell_all_conds{1}'; frac_stays_in_per_cell_all_conds{2}';frac_stays_in_per_cell_all_conds{3}';frac_stays_in_per_cell_all_conds{4}'];
no_net_frac_cell_types=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(frac_stays_in_per_cell_all_conds{1}) numel(frac_stays_in_per_cell_all_conds{2}) numel(frac_stays_in_per_cell_all_conds{3}) numel(frac_stays_in_per_cell_all_conds{4})])]';
no_net_frac_data_tbl=table(no_net_frac_values,no_net_frac_cell_types);
[no_net_frac_p,no_net_frac_tbl_anova,no_net_frac_stats]=anova1(no_net_frac_values,no_net_frac_cell_types,"on");
[no_net_frac_results,a,b,gnames]=multcompare(no_net_frac_stats);
no_net_frac_tbl = array2table(no_net_frac_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
no_net_frac_tbl.("Group A")=gnames(no_net_frac_tbl.("Group A"));
no_net_frac_tbl.("Group B")=gnames(no_net_frac_tbl.("Group B"));

tot_ent_values=[tot_enters_per_cell_all_conds{1}';tot_enters_per_cell_all_conds{2}';tot_enters_per_cell_all_conds{3}';tot_enters_per_cell_all_conds{4}'];
tot_ent_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(tot_enters_per_cell_all_conds{1}) numel(tot_enters_per_cell_all_conds{2}) numel(tot_enters_per_cell_all_conds{3}) numel(tot_enters_per_cell_all_conds{4})])]';
tot_ent_data_tbl=table(tot_ent_values,tot_ent_groups);
[tot_ent_p,tot_ent_tbl_anova,tot_ent_stats]=anova1(tot_ent_values,tot_ent_groups,"off");
[tot_ent_results,m,h,gnames]=multcompare(tot_ent_stats);
tot_ent_tbl = array2table(tot_ent_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
tot_ent_tbl.("Group A") = gnames(tot_ent_tbl.("Group A"));
tot_ent_tbl.("Group B") = gnames(tot_ent_tbl.("Group B"));

ant_ret_values=[ant_vs_ret{1}';ant_vs_ret{2}';ant_vs_ret{3}';ant_vs_ret{4}'];
ant_ret_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(ant_vs_ret{1}) numel(ant_vs_ret{2}) numel(ant_vs_ret{3}) numel(ant_vs_ret{4})])]';
ant_ret_data_tbl=table(ant_ret_values,ant_ret_groups);
[ant_ret_p,ant_ret_tbl_anova,ant_ret_stats]=anova1(ant_ret_values,ant_ret_groups,"off");
[ant_ret_results,m,h,gnames]=multcompare(ant_ret_stats);
ant_ret_tbl = array2table(ant_ret_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
ant_ret_tbl.("Group A") = gnames(ant_ret_tbl.("Group A"));
ant_ret_tbl.("Group B") = gnames(ant_ret_tbl.("Group B"));

% 
tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
FolderName_1 = tempdir_1;   % Your destination folder
% %Saving stats to a table
% writetable(frac_tbl,fullfile(FolderName_1, [fileprefix,'frac.csv']),'WriteRowNames',true);
% writetable(frac_data_tbl,fullfile(FolderName_1, [fileprefix,'frac_data.csv']),'WriteRowNames',true);
% writetable(ant_frac_tbl,fullfile(FolderName_1, [fileprefix,'ant_frac.csv']),'WriteRowNames',true);
% writetable(ant_frac_data_tbl,fullfile(FolderName_1, [fileprefix,'ant_frac_data.csv']),'WriteRowNames',true);
% writetable(ret_frac_tbl,fullfile(FolderName_1, [fileprefix,'ret_frac.csv']),'WriteRowNames',true);
% writetable(ret_frac_data_tbl,fullfile(FolderName_1, [fileprefix,'ret_frac_data.csv']),'WriteRowNames',true);
% writetable(no_net_frac_tbl,fullfile(FolderName_1, [fileprefix,'no_net_frac.csv']),'WriteRowNames',true);
% writetable(no_net_frac_data_tbl,fullfile(FolderName_1, [fileprefix,'no_net_frac_data.csv']),'WriteRowNames',true);
% writetable(tot_ent_tbl, fullfile(FolderName_1, [fileprefix,'tot_ent.csv']),'WriteRowNames',true);
% writetable(tot_ent_data_tbl, fullfile(FolderName_1, [fileprefix,'tot_ent_data.csv']),'WriteRowNames',true);
% writetable(ant_ret_tbl, fullfile(FolderName_1, [fileprefix,'ant_ret.csv']),'WriteRowNames',true);
% writetable(ant_ret_data_tbl, fullfile(FolderName_1, [fileprefix,'ant_ret_data.csv']),'WriteRowNames',true); 

FolderName = tempdir;  % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
 FigHandle = FigList(iFig);
 FigName  = get(FigHandle, 'Name');
 savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
 saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end
