set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
% polyQ length
% cmap = colormap(hot(100));
% colour_30Q=[0.5 0.5 0.5];
% colour_45Q=cmap(20,:);%[0.75 0 0];
% colour_65Q=cmap(25,:);%[0.5 0 0];
% colour_81Q=cmap(30,:);%[0.25 0 0];

%stress condition
cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_65Q=cmap(30,:);%[0.5 0 0];
cmap2= colormap(spring(100));
colour_45Q=[1.0000 0.4444 0.5556];%cmap2(50,:);%[0.75 0 0];
colour_81Q=cmap2(30,:);%[0.25 0 0];

close all;

% pth='/Volumes/Emily_htt_2/Neuron/bdnf_flux/';
pth='/Volumes/Emily_htt_2/Neuron/lyso_flux/';
% pth='/Volumes/Emily_htt_2/Neuron/mito_flux/';

fileprefix='20240119_isoHD_lyso_IFg_multi_box';

% load(fullfile(pth,'IFg_KB_frac_ant_ex_IFg_0pt25um.mat'));
% frac_ant_pt25um=frac_ant_exits_per_cell_all_conds;
% load(fullfile(pth,'IFg_KB_frac_ant_ex_IFg_0pt5um.mat'));
% frac_ant_pt5um=frac_ant_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ant_ex_IFg_1um.mat'));
frac_ant_1um=frac_ant_exits_per_cell_all_conds;
% load(fullfile(pth,'KB_frac_ant_ex_IFg_1pt5um.mat'));
% frac_ant_1pt5um=frac_ant_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ant_ex_IFg_2um.mat'));
frac_ant_2um=frac_ant_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ant_ex_IFg_5um.mat'));
frac_ant_5um=frac_ant_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ant_ex_IFg_10um.mat'));
frac_ant_10um=frac_ant_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ant_ex_IFg_15um.mat'));
frac_ant_15um=frac_ant_exits_per_cell_all_conds;

% load(fullfile(pth,'KB_frac_ret_ex_IFg_0pt25um.mat'));
% frac_ret_pt25um=frac_ret_exits_per_cell_all_conds;
% load(fullfile(pth,'KB_frac_ret_ex_IFg_0pt5um.mat'));
% frac_ret_pt5um=frac_ret_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ret_ex_IFg_1um.mat'));
frac_ret_1um=frac_ret_exits_per_cell_all_conds;
% load(fullfile(pth,'KB_frac_ret_ex_IFg_1pt5um.mat'));
% frac_ret_1pt5um=frac_ret_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ret_ex_IFg_2um.mat'));
frac_ret_2um=frac_ret_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ret_ex_IFg_5um.mat'));
frac_ret_5um=frac_ret_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ret_ex_IFg_10um.mat'));
frac_ret_10um=frac_ret_exits_per_cell_all_conds;
load(fullfile(pth,'KB_frac_ret_ex_IFg_15um.mat'));
frac_ret_15um=frac_ret_exits_per_cell_all_conds;

% load(fullfile(pth,'KB_total_enters_0pt25um.mat'));
% tot_ent_pt25um=tot_enters_per_cell_all_conds;
% load(fullfile(pth,'KB_total_enters_0pt5um.mat'));
% tot_ent_pt5um=tot_enters_per_cell_all_conds;
load(fullfile(pth,'KB_total_enters_IFg_1um.mat'));
tot_ent_1um=tot_enters_per_cell_all_conds;
% load(fullfile(pth,'KB_total_enters_1pt5um.mat'));
% tot_ent_1pt5um=tot_enters_per_cell_all_conds;
load(fullfile(pth,'KB_total_enters_IFg_2um.mat'));
tot_ent_2um=tot_enters_per_cell_all_conds;
load(fullfile(pth,'KB_total_enters_IFg_5um.mat'));
tot_ent_5um=tot_enters_per_cell_all_conds;
load(fullfile(pth,'KB_total_enters_IFg_10um.mat'));
tot_ent_10um=tot_enters_per_cell_all_conds;
load(fullfile(pth,'KB_total_enters_IFg_15um.mat'));
tot_ent_15um=tot_enters_per_cell_all_conds;

% load(fullfile(pth,'KB_frac_stay_in_0pt25um.mat'));
% frac_no_net_pt25um=frac_stays_in_per_cell_all_conds;
% load(fullfile(pth,'KB_frac_stay_in_0pt5um.mat'));
% frac_no_net_pt5um=frac_stays_in_per_cell_all_conds;
load(fullfile(pth,'KB_frac_stay_in_IFg_1um.mat'));
frac_no_net_1um=frac_stays_in_per_cell_all_conds;
% load(fullfile(pth,'KB_frac_stay_in_1pt5um.mat'));
% frac_no_net_1pt5um=frac_stays_in_per_cell_all_conds;
load(fullfile(pth,'KB_frac_stay_in_IFg_2um.mat'));
frac_no_net_2um=frac_stays_in_per_cell_all_conds;
load(fullfile(pth,'KB_frac_stay_in_IFg_5um.mat'));
frac_no_net_5um=frac_stays_in_per_cell_all_conds;
load(fullfile(pth,'KB_frac_stay_in_IFg_10um.mat'));
frac_no_net_10um=frac_stays_in_per_cell_all_conds;
load(fullfile(pth,'KB_frac_stay_in_IFg_15um.mat'));
frac_no_net_15um=frac_stays_in_per_cell_all_conds;

%Defining number of anterograde cargoes
% ant_num_pt25um{1}=frac_ant_pt25um{1}.*tot_ent_pt25um{1};
% ant_num_pt25um{2}=frac_ant_pt25um{2}.*tot_ent_pt25um{2};
% ant_num_pt25um{3}=frac_ant_pt25um{3}.*tot_ent_pt25um{3};
% ant_num_pt25um{4}=frac_ant_pt25um{4}.*tot_ent_pt25um{4};
% ant_num_pt5um{1}=frac_ant_pt5um{1}.*tot_ent_pt5um{1};
% ant_num_pt5um{2}=frac_ant_pt5um{2}.*tot_ent_pt5um{2};
% ant_num_pt5um{3}=frac_ant_pt5um{3}.*tot_ent_pt5um{3};
% ant_num_pt5um{4}=frac_ant_pt5um{4}.*tot_ent_pt5um{4};
ant_num_1um{1}=frac_ant_1um{1}.*tot_ent_1um{1};
ant_num_1um{2}=frac_ant_1um{2}.*tot_ent_1um{2};
ant_num_1um{3}=frac_ant_1um{3}.*tot_ent_1um{3};
ant_num_1um{4}=frac_ant_1um{4}.*tot_ent_1um{4};
% ant_num_1pt5um{1}=frac_ant_1pt5um{1}.*tot_ent_1pt5um{1};
% ant_num_1pt5um{2}=frac_ant_1pt5um{2}.*tot_ent_1pt5um{2};
% ant_num_1pt5um{3}=frac_ant_1pt5um{3}.*tot_ent_1pt5um{3};
% ant_num_1pt5um{4}=frac_ant_1pt5um{4}.*tot_ent_1pt5um{4};
ant_num_2um{1}=frac_ant_2um{1}.*tot_ent_2um{1};
ant_num_2um{2}=frac_ant_2um{2}.*tot_ent_2um{2};
ant_num_2um{3}=frac_ant_2um{3}.*tot_ent_2um{3};
ant_num_2um{4}=frac_ant_2um{4}.*tot_ent_2um{4};
ant_num_5um{1}=frac_ant_5um{1}.*tot_ent_5um{1};
ant_num_5um{2}=frac_ant_5um{2}.*tot_ent_5um{2};
ant_num_5um{3}=frac_ant_5um{3}.*tot_ent_5um{3};
ant_num_5um{4}=frac_ant_5um{4}.*tot_ent_5um{4};
ant_num_10um{1}=frac_ant_10um{1}.*tot_ent_10um{1};
ant_num_10um{2}=frac_ant_10um{2}.*tot_ent_10um{2};
ant_num_10um{3}=frac_ant_10um{3}.*tot_ent_10um{3};
ant_num_10um{4}=frac_ant_10um{4}.*tot_ent_10um{4};
ant_num_15um{1}=frac_ant_15um{1}.*tot_ent_15um{1};
ant_num_15um{2}=frac_ant_15um{2}.*tot_ent_15um{2};
ant_num_15um{3}=frac_ant_15um{3}.*tot_ent_15um{3};
ant_num_15um{4}=frac_ant_15um{4}.*tot_ent_15um{4};

%Defining number of retrograde cargoes
% ret_num_pt25um{1}=frac_ret_pt25um{1}.*tot_ent_pt25um{1};
% ret_num_pt25um{2}=frac_ret_pt25um{2}.*tot_ent_pt25um{2};
% ret_num_pt25um{3}=frac_ret_pt25um{3}.*tot_ent_pt25um{3};
% ret_num_pt25um{4}=frac_ret_pt25um{4}.*tot_ent_pt25um{4};
% ret_num_pt5um{1}=frac_ret_pt5um{1}.*tot_ent_pt5um{1};
% ret_num_pt5um{2}=frac_ret_pt5um{2}.*tot_ent_pt5um{2};
% ret_num_pt5um{3}=frac_ret_pt5um{3}.*tot_ent_pt5um{3};
% ret_num_pt5um{4}=frac_ret_pt5um{4}.*tot_ent_pt5um{4};
ret_num_1um{1}=frac_ret_1um{1}.*tot_ent_1um{1};
ret_num_1um{2}=frac_ret_1um{2}.*tot_ent_1um{2};
ret_num_1um{3}=frac_ret_1um{3}.*tot_ent_1um{3};
ret_num_1um{4}=frac_ret_1um{4}.*tot_ent_1um{4};
% ret_num_1pt5um{1}=frac_ret_1pt5um{1}.*tot_ent_1pt5um{1};
% ret_num_1pt5um{2}=frac_ret_1pt5um{2}.*tot_ent_1pt5um{2};
% ret_num_1pt5um{3}=frac_ret_1pt5um{3}.*tot_ent_1pt5um{3};
% ret_num_1pt5um{4}=frac_ret_1pt5um{4}.*tot_ent_1pt5um{4};
ret_num_2um{1}=frac_ret_2um{1}.*tot_ent_2um{1};
ret_num_2um{2}=frac_ret_2um{2}.*tot_ent_2um{2};
ret_num_2um{3}=frac_ret_2um{3}.*tot_ent_2um{3};
ret_num_2um{4}=frac_ret_2um{4}.*tot_ent_2um{4};
ret_num_5um{1}=frac_ret_5um{1}.*tot_ent_5um{1};
ret_num_5um{2}=frac_ret_5um{2}.*tot_ent_5um{2};
ret_num_5um{3}=frac_ret_5um{3}.*tot_ent_5um{3};
ret_num_5um{4}=frac_ret_5um{4}.*tot_ent_5um{4};
ret_num_10um{1}=frac_ret_10um{1}.*tot_ent_10um{1};
ret_num_10um{2}=frac_ret_10um{2}.*tot_ent_10um{2};
ret_num_10um{3}=frac_ret_10um{3}.*tot_ent_10um{3};
ret_num_10um{4}=frac_ret_10um{4}.*tot_ent_10um{4};
ret_num_15um{1}=frac_ret_15um{1}.*tot_ent_15um{1};
ret_num_15um{2}=frac_ret_15um{2}.*tot_ent_15um{2};
ret_num_15um{3}=frac_ret_15um{3}.*tot_ent_15um{3};
ret_num_15um{4}=frac_ret_15um{4}.*tot_ent_15um{4};

%Defining number of cargoes with no net flux
% no_net_num_pt25um{1}=frac_no_net_pt25um{1}.*tot_ent_pt25um{1};
% no_net_num_pt25um{2}=frac_no_net_pt25um{2}.*tot_ent_pt25um{2};
% no_net_num_pt25um{3}=frac_no_net_pt25um{3}.*tot_ent_pt25um{3};
% no_net_num_pt25um{4}=frac_no_net_pt25um{4}.*tot_ent_pt25um{4};
% no_net_num_pt5um{1}=frac_no_net_pt5um{1}.*tot_ent_pt5um{1};
% no_net_num_pt5um{2}=frac_no_net_pt5um{2}.*tot_ent_pt5um{2};
% no_net_num_pt5um{3}=frac_no_net_pt5um{3}.*tot_ent_pt5um{3};
% no_net_num_pt5um{4}=frac_no_net_pt5um{4}.*tot_ent_pt5um{4};
no_net_num_1um{1}=frac_no_net_1um{1}.*tot_ent_1um{1};
no_net_num_1um{2}=frac_no_net_1um{2}.*tot_ent_1um{2};
no_net_num_1um{3}=frac_no_net_1um{3}.*tot_ent_1um{3};
no_net_num_1um{4}=frac_no_net_1um{4}.*tot_ent_1um{4};
% no_net_num_1pt5um{1}=frac_no_net_1pt5um{1}.*tot_ent_1pt5um{1};
% no_net_num_1pt5um{2}=frac_no_net_1pt5um{2}.*tot_ent_1pt5um{2};
% no_net_num_1pt5um{3}=frac_no_net_1pt5um{3}.*tot_ent_1pt5um{3};
% no_net_num_1pt5um{4}=frac_no_net_1pt5um{4}.*tot_ent_1pt5um{4};
no_net_num_2um{1}=frac_no_net_2um{1}.*tot_ent_2um{1};
no_net_num_2um{2}=frac_no_net_2um{2}.*tot_ent_2um{2};
no_net_num_2um{3}=frac_no_net_2um{3}.*tot_ent_2um{3};
no_net_num_2um{4}=frac_no_net_2um{4}.*tot_ent_2um{4};
no_net_num_5um{1}=frac_no_net_5um{1}.*tot_ent_5um{1};
no_net_num_5um{2}=frac_no_net_5um{2}.*tot_ent_5um{2};
no_net_num_5um{3}=frac_no_net_5um{3}.*tot_ent_5um{3};
no_net_num_5um{4}=frac_no_net_5um{4}.*tot_ent_5um{4};
no_net_num_10um{1}=frac_no_net_10um{1}.*tot_ent_10um{1};
no_net_num_10um{2}=frac_no_net_10um{2}.*tot_ent_10um{2};
no_net_num_10um{3}=frac_no_net_10um{3}.*tot_ent_10um{3};
no_net_num_10um{4}=frac_no_net_10um{4}.*tot_ent_10um{4};
no_net_num_15um{1}=frac_no_net_15um{1}.*tot_ent_15um{1};
no_net_num_15um{2}=frac_no_net_15um{2}.*tot_ent_15um{2};
no_net_num_15um{3}=frac_no_net_15um{3}.*tot_ent_15um{3};
no_net_num_15um{4}=frac_no_net_15um{4}.*tot_ent_15um{4};

%Total enters means
lb_mean_30Q_tot_enters=[mean(tot_ent_1um{1}),mean(tot_ent_2um{1}),mean(tot_ent_5um{1}),mean(tot_ent_10um{1}),mean(tot_ent_15um{1})];
lb_mean_45Q_tot_enters=[mean(tot_ent_1um{2}),mean(tot_ent_2um{2}),mean(tot_ent_5um{2}),mean(tot_ent_10um{2}),mean(tot_ent_15um{2})];
lb_mean_65Q_tot_enters=[mean(tot_ent_1um{3}),mean(tot_ent_2um{3}),mean(tot_ent_5um{3}),mean(tot_ent_10um{3}),mean(tot_ent_15um{3})];
lb_mean_81Q_tot_enters=[mean(tot_ent_1um{4}),mean(tot_ent_2um{4}),mean(tot_ent_5um{4}),mean(tot_ent_10um{4}),mean(tot_ent_15um{4})];

lb_mean_30Q_tot_enters_pm=[mean(tot_ent_1um{1}./5),mean(tot_ent_2um{1}./10),mean(tot_ent_5um{1}./25),mean(tot_ent_10um{1}./50),mean(tot_ent_15um{1}./75)];
lb_mean_45Q_tot_enters_pm=[mean(tot_ent_1um{2}./5),mean(tot_ent_2um{2}./10),mean(tot_ent_5um{2}./25),mean(tot_ent_10um{2}./50),mean(tot_ent_15um{2}./75)];
lb_mean_65Q_tot_enters_pm=[mean(tot_ent_1um{3}./5),mean(tot_ent_2um{3}./10),mean(tot_ent_5um{3}./25),mean(tot_ent_10um{3}./50),mean(tot_ent_15um{3}./75)];
lb_mean_81Q_tot_enters_pm=[mean(tot_ent_1um{4}./5),mean(tot_ent_2um{4}./10),mean(tot_ent_5um{4}./25),mean(tot_ent_10um{4}./50),mean(tot_ent_15um{4}./75)];


% all_mean_30Q_tot_enters=[mean(tot_ent_pt25um{1},'omitnan'),mean(tot_ent_pt5um{1},'omitnan'),mean(tot_ent_1um{1},'omitnan'),mean(tot_ent_1pt5um{1},'omitnan'),mean(tot_ent_2um{1},'omitnan'),mean(tot_ent_5um{1},'omitnan'),mean(tot_ent_10um{1},'omitnan'),mean(tot_ent_15um{1},'omitnan')];
% all_mean_45Q_tot_enters=[mean(tot_ent_pt25um{2},'omitnan'),mean(tot_ent_pt5um{2},'omitnan'),mean(tot_ent_1um{2},'omitnan'),mean(tot_ent_1pt5um{2},'omitnan'),mean(tot_ent_2um{2},'omitnan'),mean(tot_ent_5um{2},'omitnan'),mean(tot_ent_10um{2},'omitnan'),mean(tot_ent_15um{2},'omitnan')];
% all_mean_65Q_tot_enters=[mean(tot_ent_pt25um{3},'omitnan'),mean(tot_ent_pt5um{3},'omitnan'),mean(tot_ent_1um{3},'omitnan'),mean(tot_ent_1pt5um{3},'omitnan'),mean(tot_ent_2um{3},'omitnan'),mean(tot_ent_5um{3},'omitnan'),mean(tot_ent_10um{3},'omitnan'),mean(tot_ent_15um{3},'omitnan')];
% all_mean_81Q_tot_enters=[mean(tot_ent_pt25um{4},'omitnan'),mean(tot_ent_pt5um{4},'omitnan'),mean(tot_ent_1um{4},'omitnan'),mean(tot_ent_1pt5um{4},'omitnan'),mean(tot_ent_2um{4},'omitnan'),mean(tot_ent_5um{4},'omitnan'),mean(tot_ent_10um{4},'omitnan'),mean(tot_ent_15um{4},'omitnan')];

%Total enters: calculating SEM for errorbars
lb_SEM_30Q_tot_enters{1}=std(tot_ent_1um{1},'omitnan')./(sqrt(length(tot_ent_1um{1})));
lb_SEM_30Q_tot_enters{2}=std(tot_ent_2um{1},'omitnan')./(sqrt(length(tot_ent_2um{1})));
lb_SEM_30Q_tot_enters{3}=std(tot_ent_5um{1},'omitnan')./(sqrt(length(tot_ent_5um{1})));
lb_SEM_30Q_tot_enters{4}=std(tot_ent_10um{1},'omitnan')./(sqrt(length(tot_ent_10um{1})));
lb_SEM_30Q_tot_enters{5}=std(tot_ent_15um{1},'omitnan')./(sqrt(length(tot_ent_15um{1})));

lb_SEM_30Q_tot_enters_pm{1}=std(tot_ent_1um{1}./5,'omitnan')./(sqrt(length(tot_ent_1um{1})));
lb_SEM_30Q_tot_enters_pm{2}=std(tot_ent_2um{1}./10,'omitnan')./(sqrt(length(tot_ent_2um{1})));
lb_SEM_30Q_tot_enters_pm{3}=std(tot_ent_5um{1}./25,'omitnan')./(sqrt(length(tot_ent_5um{1})));
lb_SEM_30Q_tot_enters_pm{4}=std(tot_ent_10um{1}./50,'omitnan')./(sqrt(length(tot_ent_10um{1})));
lb_SEM_30Q_tot_enters_pm{5}=std(tot_ent_15um{1}./75,'omitnan')./(sqrt(length(tot_ent_15um{1})));

% all_SEM_30Q_tot_enters{1}=std(tot_ent_pt25um{1},'omitnan')./(sqrt(length(tot_ent_pt25um{1})));
% all_SEM_30Q_tot_enters{2}=std(tot_ent_pt5um{1},'omitnan')./(sqrt(length(tot_ent_pt5um{1})));
% all_SEM_30Q_tot_enters{3}=std(tot_ent_1um{1},'omitnan')./(sqrt(length(tot_ent_1um{1})));
% all_SEM_30Q_tot_enters{4}=std(tot_ent_1pt5um{1},'omitnan')./(sqrt(length(tot_ent_1pt5um{1})));
% all_SEM_30Q_tot_enters{5}=std(tot_ent_2um{1},'omitnan')./(sqrt(length(tot_ent_2um{1})));
% all_SEM_30Q_tot_enters{6}=std(tot_ent_5um{1},'omitnan')./(sqrt(length(tot_ent_5um{1})));
% all_SEM_30Q_tot_enters{7}=std(tot_ent_10um{1},'omitnan')./(sqrt(length(tot_ent_10um{1})));
% all_SEM_30Q_tot_enters{8}=std(tot_ent_15um{1},'omitnan')./(sqrt(length(tot_ent_15um{1})));

lb_SEM_45Q_tot_enters{1}=std(tot_ent_1um{2},'omitnan')./(sqrt(length(tot_ent_1um{2})));
lb_SEM_45Q_tot_enters{2}=std(tot_ent_2um{2},'omitnan')./(sqrt(length(tot_ent_2um{2})));
lb_SEM_45Q_tot_enters{3}=std(tot_ent_5um{2},'omitnan')./(sqrt(length(tot_ent_5um{2})));
lb_SEM_45Q_tot_enters{4}=std(tot_ent_10um{2},'omitnan')./(sqrt(length(tot_ent_10um{2})));
lb_SEM_45Q_tot_enters{5}=std(tot_ent_15um{2},'omitnan')./(sqrt(length(tot_ent_15um{2})));

lb_SEM_45Q_tot_enters_pm{1}=std(tot_ent_1um{2}./5,'omitnan')./(sqrt(length(tot_ent_1um{2})));
lb_SEM_45Q_tot_enters_pm{2}=std(tot_ent_2um{2}./10,'omitnan')./(sqrt(length(tot_ent_2um{2})));
lb_SEM_45Q_tot_enters_pm{3}=std(tot_ent_5um{2}./25,'omitnan')./(sqrt(length(tot_ent_5um{2})));
lb_SEM_45Q_tot_enters_pm{4}=std(tot_ent_10um{2}./50,'omitnan')./(sqrt(length(tot_ent_10um{2})));
lb_SEM_45Q_tot_enters_pm{5}=std(tot_ent_15um{2}./75,'omitnan')./(sqrt(length(tot_ent_15um{2})));

% all_SEM_45Q_tot_enters{1}=std(tot_ent_pt25um{2},'omitnan')./(sqrt(length(tot_ent_pt25um{2})));
% all_SEM_45Q_tot_enters{2}=std(tot_ent_pt5um{2},'omitnan')./(sqrt(length(tot_ent_pt5um{2})));
% all_SEM_45Q_tot_enters{3}=std(tot_ent_1um{2},'omitnan')./(sqrt(length(tot_ent_1um{2})));
% all_SEM_45Q_tot_enters{4}=std(tot_ent_1pt5um{2},'omitnan')./(sqrt(length(tot_ent_1pt5um{2})));
% all_SEM_45Q_tot_enters{5}=std(tot_ent_2um{2},'omitnan')./(sqrt(length(tot_ent_2um{2})));
% all_SEM_45Q_tot_enters{6}=std(tot_ent_5um{2},'omitnan')./(sqrt(length(tot_ent_5um{2})));
% all_SEM_45Q_tot_enters{7}=std(tot_ent_10um{2},'omitnan')./(sqrt(length(tot_ent_10um{2})));
% all_SEM_45Q_tot_enters{8}=std(tot_ent_15um{2},'omitnan')./(sqrt(length(tot_ent_15um{2})));

lb_SEM_65Q_tot_enters{1}=std(tot_ent_1um{3},'omitnan')./(sqrt(length(tot_ent_1um{3})));
lb_SEM_65Q_tot_enters{2}=std(tot_ent_2um{3},'omitnan')./(sqrt(length(tot_ent_2um{3})));
lb_SEM_65Q_tot_enters{3}=std(tot_ent_5um{3},'omitnan')./(sqrt(length(tot_ent_5um{3})));
lb_SEM_65Q_tot_enters{4}=std(tot_ent_10um{3},'omitnan')./(sqrt(length(tot_ent_10um{3})));
lb_SEM_65Q_tot_enters{5}=std(tot_ent_15um{3},'omitnan')./(sqrt(length(tot_ent_15um{3})));

lb_SEM_65Q_tot_enters_pm{1}=std(tot_ent_1um{3}./5,'omitnan')./(sqrt(length(tot_ent_1um{3})));
lb_SEM_65Q_tot_enters_pm{2}=std(tot_ent_2um{3}./10,'omitnan')./(sqrt(length(tot_ent_2um{3})));
lb_SEM_65Q_tot_enters_pm{3}=std(tot_ent_5um{3}./25,'omitnan')./(sqrt(length(tot_ent_5um{3})));
lb_SEM_65Q_tot_enters_pm{4}=std(tot_ent_10um{3}./50,'omitnan')./(sqrt(length(tot_ent_10um{3})));
lb_SEM_65Q_tot_enters_pm{5}=std(tot_ent_15um{3}./75,'omitnan')./(sqrt(length(tot_ent_15um{3})));

% all_SEM_65Q_tot_enters{1}=std(tot_ent_pt25um{3},'omitnan')./(sqrt(length(tot_ent_pt25um{3})));
% all_SEM_65Q_tot_enters{2}=std(tot_ent_pt5um{3},'omitnan')./(sqrt(length(tot_ent_pt5um{3})));
% all_SEM_65Q_tot_enters{3}=std(tot_ent_1um{3},'omitnan')./(sqrt(length(tot_ent_1um{3})));
% all_SEM_65Q_tot_enters{4}=std(tot_ent_1pt5um{3},'omitnan')./(sqrt(length(tot_ent_1pt5um{3})));
% all_SEM_65Q_tot_enters{5}=std(tot_ent_2um{3},'omitnan')./(sqrt(length(tot_ent_2um{3})));
% all_SEM_65Q_tot_enters{6}=std(tot_ent_5um{3},'omitnan')./(sqrt(length(tot_ent_5um{3})));
% all_SEM_65Q_tot_enters{7}=std(tot_ent_10um{3},'omitnan')./(sqrt(length(tot_ent_10um{3})));
% all_SEM_65Q_tot_enters{8}=std(tot_ent_15um{3},'omitnan')./(sqrt(length(tot_ent_15um{3})));

lb_SEM_81Q_tot_enters{1}=std(tot_ent_1um{4},'omitnan')./(sqrt(length(tot_ent_1um{4})));
lb_SEM_81Q_tot_enters{2}=std(tot_ent_2um{4},'omitnan')./(sqrt(length(tot_ent_2um{4})));
lb_SEM_81Q_tot_enters{3}=std(tot_ent_5um{4},'omitnan')./(sqrt(length(tot_ent_5um{4})));
lb_SEM_81Q_tot_enters{4}=std(tot_ent_10um{4},'omitnan')./(sqrt(length(tot_ent_10um{4})));
lb_SEM_81Q_tot_enters{5}=std(tot_ent_15um{4},'omitnan')./(sqrt(length(tot_ent_15um{4})));

lb_SEM_81Q_tot_enters_pm{1}=std(tot_ent_1um{4}./5,'omitnan')./(sqrt(length(tot_ent_1um{4})));
lb_SEM_81Q_tot_enters_pm{2}=std(tot_ent_2um{4}./10,'omitnan')./(sqrt(length(tot_ent_2um{4})));
lb_SEM_81Q_tot_enters_pm{3}=std(tot_ent_5um{4}./25,'omitnan')./(sqrt(length(tot_ent_5um{4})));
lb_SEM_81Q_tot_enters_pm{4}=std(tot_ent_10um{4}./50,'omitnan')./(sqrt(length(tot_ent_10um{4})));
lb_SEM_81Q_tot_enters_pm{5}=std(tot_ent_15um{4}./75,'omitnan')./(sqrt(length(tot_ent_15um{4})));

% all_SEM_81Q_tot_enters{1}=std(tot_ent_pt25um{4},'omitnan')./(sqrt(length(tot_ent_pt25um{4})));
% all_SEM_81Q_tot_enters{2}=std(tot_ent_pt5um{4},'omitnan')./(sqrt(length(tot_ent_pt5um{4})));
% all_SEM_81Q_tot_enters{3}=std(tot_ent_1um{4},'omitnan')./(sqrt(length(tot_ent_1um{4})));
% all_SEM_81Q_tot_enters{4}=std(tot_ent_1pt5um{4},'omitnan')./(sqrt(length(tot_ent_1pt5um{4})));
% all_SEM_81Q_tot_enters{5}=std(tot_ent_2um{4},'omitnan')./(sqrt(length(tot_ent_2um{4})));
% all_SEM_81Q_tot_enters{6}=std(tot_ent_5um{4},'omitnan')./(sqrt(length(tot_ent_5um{4})));
% all_SEM_81Q_tot_enters{7}=std(tot_ent_10um{4},'omitnan')./(sqrt(length(tot_ent_10um{4})));
% all_SEM_81Q_tot_enters{8}=std(tot_ent_15um{4},'omitnan')./(sqrt(length(tot_ent_15um{4})));

% %Anterograde Fraction Means
% all_mean_30Q_ant_frac=[mean(frac_ant_1um{1},'omitnan'), mean(frac_ant_2um{1},'omitnan'),mean(frac_ant_5um{1},'omitnan'),mean(frac_ant_10um{1},'omitnan'),mean(frac_ant_15um{1},'omitnan')];
% all_mean_45Q_ant_frac=[mean(frac_ant_1um{2},'omitnan'), mean(frac_ant_2um{2},'omitnan'),mean(frac_ant_5um{2},'omitnan'),mean(frac_ant_10um{2},'omitnan'),mean(frac_ant_15um{2},'omitnan')];
% all_mean_65Q_ant_frac=[mean(frac_ant_1um{3},'omitnan'), mean(frac_ant_2um{3},'omitnan'),mean(frac_ant_5um{3},'omitnan'),mean(frac_ant_10um{3},'omitnan'),mean(frac_ant_15um{3},'omitnan')];
% all_mean_81Q_ant_frac=[mean(frac_ant_1um{4},'omitnan'), mean(frac_ant_2um{4},'omitnan'),mean(frac_ant_5um{4},'omitnan'),mean(frac_ant_10um{4},'omitnan'),mean(frac_ant_15um{4},'omitnan')];
% 
% %Anterograde fraction: calculating SEM for errorbars
% all_SEM_30Q_ant_frac{1}=std(frac_ant_1um{1})./(sqrt(length(frac_ant_1um{1})));
% all_SEM_30Q_ant_frac{2}=std(frac_ant_2um{1})./(sqrt(length(frac_ant_2um{1})));
% all_SEM_30Q_ant_frac{3}=std(frac_ant_5um{1})./(sqrt(length(frac_ant_5um{1})));
% all_SEM_30Q_ant_frac{4}=std(frac_ant_10um{1})./(sqrt(length(frac_ant_10um{1})));
% all_SEM_30Q_ant_frac{5}=std(frac_ant_15um{1})./(sqrt(length(frac_ant_15um{1})));
% 
% all_SEM_45Q_ant_frac{1}=std(frac_ant_1um{2})./(sqrt(length(frac_ant_1um{2})));
% all_SEM_45Q_ant_frac{2}=std(frac_ant_2um{2})./(sqrt(length(frac_ant_2um{2})));
% all_SEM_45Q_ant_frac{3}=std(frac_ant_5um{2})./(sqrt(length(frac_ant_5um{2})));
% all_SEM_45Q_ant_frac{4}=std(frac_ant_10um{2})./(sqrt(length(frac_ant_10um{2})));
% all_SEM_45Q_ant_frac{5}=std(frac_ant_15um{2})./(sqrt(length(frac_ant_15um{2})));
% 
% all_SEM_65Q_ant_frac{1}=std(frac_ant_1um{3})./(sqrt(length(frac_ant_1um{3})));
% all_SEM_65Q_ant_frac{2}=std(frac_ant_2um{3})./(sqrt(length(frac_ant_2um{3})));
% all_SEM_65Q_ant_frac{3}=std(frac_ant_5um{3})./(sqrt(length(frac_ant_5um{3})));
% all_SEM_65Q_ant_frac{4}=std(frac_ant_10um{3})./(sqrt(length(frac_ant_10um{3})));
% all_SEM_65Q_ant_frac{5}=std(frac_ant_15um{3})./(sqrt(length(frac_ant_15um{3})));
% 
% all_SEM_81Q_ant_frac{1}=std(frac_ant_1um{4})./(sqrt(length(frac_ant_1um{4})));
% all_SEM_81Q_ant_frac{2}=std(frac_ant_2um{4})./(sqrt(length(frac_ant_2um{4})));
% all_SEM_81Q_ant_frac{3}=std(frac_ant_5um{4})./(sqrt(length(frac_ant_5um{4})));
% all_SEM_81Q_ant_frac{4}=std(frac_ant_10um{4})./(sqrt(length(frac_ant_10um{4})));
% all_SEM_81Q_ant_frac{5}=std(frac_ant_15um{4})./(sqrt(length(frac_ant_15um{4})));

%Anterograde Mean Number of Cargoes
lb_mean_30Q_ant_num=[mean(ant_num_1um{1},'omitnan'), mean(ant_num_2um{1},'omitnan'),mean(ant_num_5um{1},'omitnan'),mean(ant_num_10um{1},'omitnan'),mean(ant_num_15um{1},'omitnan')];
lb_mean_45Q_ant_num=[mean(ant_num_1um{2},'omitnan'), mean(ant_num_2um{2},'omitnan'),mean(ant_num_5um{2},'omitnan'),mean(ant_num_10um{2},'omitnan'),mean(ant_num_15um{2},'omitnan')];
lb_mean_65Q_ant_num=[mean(ant_num_1um{3},'omitnan'), mean(ant_num_2um{3},'omitnan'),mean(ant_num_5um{3},'omitnan'),mean(ant_num_10um{3},'omitnan'),mean(ant_num_15um{3},'omitnan')];
lb_mean_81Q_ant_num=[mean(ant_num_1um{4},'omitnan'), mean(ant_num_2um{4},'omitnan'),mean(ant_num_5um{4},'omitnan'),mean(ant_num_10um{4},'omitnan'),mean(ant_num_15um{4},'omitnan')];

lb_mean_30Q_ant_num_pm=[mean(ant_num_1um{1}./5,'omitnan'), mean(ant_num_2um{1}./10,'omitnan'),mean(ant_num_5um{1}./25,'omitnan'),mean(ant_num_10um{1}./50,'omitnan'),mean(ant_num_15um{1}./75,'omitnan')];
lb_mean_45Q_ant_num_pm=[mean(ant_num_1um{2}./5,'omitnan'), mean(ant_num_2um{2}./10,'omitnan'),mean(ant_num_5um{2}./25,'omitnan'),mean(ant_num_10um{2}./50,'omitnan'),mean(ant_num_15um{2}./75,'omitnan')];
lb_mean_65Q_ant_num_pm=[mean(ant_num_1um{3}./5,'omitnan'), mean(ant_num_2um{3}./10,'omitnan'),mean(ant_num_5um{3}./25,'omitnan'),mean(ant_num_10um{3}./50,'omitnan'),mean(ant_num_15um{3}./75,'omitnan')];
lb_mean_81Q_ant_num_pm=[mean(ant_num_1um{4}./5,'omitnan'), mean(ant_num_2um{4}./10,'omitnan'),mean(ant_num_5um{4}./25,'omitnan'),mean(ant_num_10um{4}./50,'omitnan'),mean(ant_num_15um{4}./75,'omitnan')];


% all_mean_30Q_ant_num=[mean(ant_num_pt25um{1},'omitnan'),mean(ant_num_pt5um{1},'omitnan'),mean(ant_num_1um{1},'omitnan'),mean(ant_num_1pt5um{1},'omitnan'),mean(ant_num_2um{1},'omitnan'),mean(ant_num_5um{1},'omitnan'),mean(ant_num_10um{1},'omitnan'),mean(ant_num_15um{1},'omitnan')];
% all_mean_45Q_ant_num=[mean(ant_num_pt25um{2},'omitnan'),mean(ant_num_pt5um{2},'omitnan'),mean(ant_num_1um{2},'omitnan'),mean(ant_num_1pt5um{2},'omitnan'),mean(ant_num_2um{2},'omitnan'),mean(ant_num_5um{2},'omitnan'),mean(ant_num_10um{2},'omitnan'),mean(ant_num_15um{2},'omitnan')];
% all_mean_65Q_ant_num=[mean(ant_num_pt25um{3},'omitnan'),mean(ant_num_pt5um{3},'omitnan'),mean(ant_num_1um{3},'omitnan'),mean(ant_num_1pt5um{3},'omitnan'),mean(ant_num_2um{3},'omitnan'),mean(ant_num_5um{3},'omitnan'),mean(ant_num_10um{3},'omitnan'),mean(ant_num_15um{3},'omitnan')];
% all_mean_81Q_ant_num=[mean(ant_num_pt25um{4},'omitnan'),mean(ant_num_pt5um{4},'omitnan'),mean(ant_num_1um{4},'omitnan'),mean(ant_num_1pt5um{4},'omitnan'),mean(ant_num_2um{4},'omitnan'),mean(ant_num_5um{4},'omitnan'),mean(ant_num_10um{4},'omitnan'),mean(ant_num_15um{4},'omitnan')];

%Anterograde number: calculating SEM for errorbars
lb_SEM_30Q_ant_num{1}=std(ant_num_1um{1},'omitnan')./(sqrt(length(ant_num_1um{1})));
lb_SEM_30Q_ant_num{2}=std(ant_num_2um{1},'omitnan')./(sqrt(length(ant_num_2um{1})));
lb_SEM_30Q_ant_num{3}=std(ant_num_5um{1},'omitnan')./(sqrt(length(ant_num_5um{1})));
lb_SEM_30Q_ant_num{4}=std(ant_num_10um{1},'omitnan')./(sqrt(length(ant_num_10um{1})));
lb_SEM_30Q_ant_num{5}=std(ant_num_15um{1},'omitnan')./(sqrt(length(ant_num_15um{1})));

lb_SEM_30Q_ant_num_pm{1}=std(ant_num_1um{1}./5,'omitnan')./(sqrt(length(ant_num_1um{1})));
lb_SEM_30Q_ant_num_pm{2}=std(ant_num_2um{1}./10,'omitnan')./(sqrt(length(ant_num_2um{1})));
lb_SEM_30Q_ant_num_pm{3}=std(ant_num_5um{1}./25,'omitnan')./(sqrt(length(ant_num_5um{1})));
lb_SEM_30Q_ant_num_pm{4}=std(ant_num_10um{1}./50,'omitnan')./(sqrt(length(ant_num_10um{1})));
lb_SEM_30Q_ant_num_pm{5}=std(ant_num_15um{1}./75,'omitnan')./(sqrt(length(ant_num_15um{1})));

% all_SEM_30Q_ant_num{1}=std(ant_num_pt25um{1},'omitnan')./(sqrt(length(ant_num_pt25um{1})));
% all_SEM_30Q_ant_num{2}=std(ant_num_pt5um{1},'omitnan')./(sqrt(length(ant_num_pt5um{1})));
% all_SEM_30Q_ant_num{3}=std(ant_num_1um{1},'omitnan')./(sqrt(length(ant_num_1um{1})));
% all_SEM_30Q_ant_num{4}=std(ant_num_1pt5um{1},'omitnan')./(sqrt(length(ant_num_1pt5um{1})));
% all_SEM_30Q_ant_num{5}=std(ant_num_2um{1},'omitnan')./(sqrt(length(ant_num_2um{1})));
% all_SEM_30Q_ant_num{6}=std(ant_num_5um{1},'omitnan')./(sqrt(length(ant_num_5um{1})));
% all_SEM_30Q_ant_num{7}=std(ant_num_10um{1},'omitnan')./(sqrt(length(ant_num_10um{1})));
% all_SEM_30Q_ant_num{8}=std(ant_num_15um{1},'omitnan')./(sqrt(length(ant_num_15um{1})));

lb_SEM_45Q_ant_num{1}=std(ant_num_1um{2},'omitnan')./(sqrt(length(ant_num_1um{2})));
lb_SEM_45Q_ant_num{2}=std(ant_num_2um{2},'omitnan')./(sqrt(length(ant_num_2um{2})));
lb_SEM_45Q_ant_num{3}=std(ant_num_5um{2},'omitnan')./(sqrt(length(ant_num_5um{2})));
lb_SEM_45Q_ant_num{4}=std(ant_num_10um{2},'omitnan')./(sqrt(length(ant_num_10um{2})));
lb_SEM_45Q_ant_num{5}=std(ant_num_15um{2},'omitnan')./(sqrt(length(ant_num_15um{2})));

lb_SEM_45Q_ant_num_pm{1}=std(ant_num_1um{2}./5,'omitnan')./(sqrt(length(ant_num_1um{2})));
lb_SEM_45Q_ant_num_pm{2}=std(ant_num_2um{2}./10,'omitnan')./(sqrt(length(ant_num_2um{2})));
lb_SEM_45Q_ant_num_pm{3}=std(ant_num_5um{2}./25,'omitnan')./(sqrt(length(ant_num_5um{2})));
lb_SEM_45Q_ant_num_pm{4}=std(ant_num_10um{2}./50,'omitnan')./(sqrt(length(ant_num_10um{2})));
lb_SEM_45Q_ant_num_pm{5}=std(ant_num_15um{2}./75,'omitnan')./(sqrt(length(ant_num_15um{2})));

% all_SEM_45Q_ant_num{1}=std(ant_num_pt25um{2},'omitnan')./(sqrt(length(ant_num_pt25um{2})));
% all_SEM_45Q_ant_num{2}=std(ant_num_pt5um{2},'omitnan')./(sqrt(length(ant_num_pt5um{2})));
% all_SEM_45Q_ant_num{3}=std(ant_num_1um{2},'omitnan')./(sqrt(length(ant_num_1um{2})));
% all_SEM_45Q_ant_num{4}=std(ant_num_1pt5um{2},'omitnan')./(sqrt(length(ant_num_1pt5um{2})));
% all_SEM_45Q_ant_num{5}=std(ant_num_2um{2},'omitnan')./(sqrt(length(ant_num_2um{2})));
% all_SEM_45Q_ant_num{6}=std(ant_num_5um{2},'omitnan')./(sqrt(length(ant_num_5um{2})));
% all_SEM_45Q_ant_num{7}=std(ant_num_10um{2},'omitnan')./(sqrt(length(ant_num_10um{2})));
% all_SEM_45Q_ant_num{8}=std(ant_num_15um{2},'omitnan')./(sqrt(length(ant_num_15um{2})));

lb_SEM_65Q_ant_num{1}=std(ant_num_1um{3},'omitnan')./(sqrt(length(ant_num_1um{3})));
lb_SEM_65Q_ant_num{2}=std(ant_num_2um{3},'omitnan')./(sqrt(length(ant_num_2um{3})));
lb_SEM_65Q_ant_num{3}=std(ant_num_5um{3},'omitnan')./(sqrt(length(ant_num_5um{3})));
lb_SEM_65Q_ant_num{4}=std(ant_num_10um{3},'omitnan')./(sqrt(length(ant_num_10um{3})));
lb_SEM_65Q_ant_num{5}=std(ant_num_15um{3},'omitnan')./(sqrt(length(ant_num_15um{3})));

lb_SEM_65Q_ant_num_pm{1}=std(ant_num_1um{3}./5,'omitnan')./(sqrt(length(ant_num_1um{3})));
lb_SEM_65Q_ant_num_pm{2}=std(ant_num_2um{3}./10,'omitnan')./(sqrt(length(ant_num_2um{3})));
lb_SEM_65Q_ant_num_pm{3}=std(ant_num_5um{3}./25,'omitnan')./(sqrt(length(ant_num_5um{3})));
lb_SEM_65Q_ant_num_pm{4}=std(ant_num_10um{3}./50,'omitnan')./(sqrt(length(ant_num_10um{3})));
lb_SEM_65Q_ant_num_pm{5}=std(ant_num_15um{3}./75,'omitnan')./(sqrt(length(ant_num_15um{3})));

% all_SEM_65Q_ant_num{1}=std(ant_num_pt25um{3},'omitnan')./(sqrt(length(ant_num_pt25um{3})));
% all_SEM_65Q_ant_num{2}=std(ant_num_pt5um{3},'omitnan')./(sqrt(length(ant_num_pt5um{3})));
% all_SEM_65Q_ant_num{3}=std(ant_num_1um{3},'omitnan')./(sqrt(length(ant_num_1um{3})));
% all_SEM_65Q_ant_num{4}=std(ant_num_1pt5um{3},'omitnan')./(sqrt(length(ant_num_1pt5um{3})));
% all_SEM_65Q_ant_num{5}=std(ant_num_2um{3},'omitnan')./(sqrt(length(ant_num_2um{3})));
% all_SEM_65Q_ant_num{6}=std(ant_num_5um{3},'omitnan')./(sqrt(length(ant_num_5um{3})));
% all_SEM_65Q_ant_num{7}=std(ant_num_10um{3},'omitnan')./(sqrt(length(ant_num_10um{3})));
% all_SEM_65Q_ant_num{8}=std(ant_num_15um{3},'omitnan')./(sqrt(length(ant_num_15um{3})));

lb_SEM_81Q_ant_num{1}=std(ant_num_1um{4},'omitnan')./(sqrt(length(ant_num_1um{4})));
lb_SEM_81Q_ant_num{2}=std(ant_num_2um{4},'omitnan')./(sqrt(length(ant_num_2um{4})));
lb_SEM_81Q_ant_num{3}=std(ant_num_5um{4},'omitnan')./(sqrt(length(ant_num_5um{4})));
lb_SEM_81Q_ant_num{4}=std(ant_num_10um{4},'omitnan')./(sqrt(length(ant_num_10um{4})));
lb_SEM_81Q_ant_num{5}=std(ant_num_15um{4},'omitnan')./(sqrt(length(ant_num_15um{4})));

lb_SEM_81Q_ant_num_pm{1}=std(ant_num_1um{4}./5,'omitnan')./(sqrt(length(ant_num_1um{4})));
lb_SEM_81Q_ant_num_pm{2}=std(ant_num_2um{4}./10,'omitnan')./(sqrt(length(ant_num_2um{4})));
lb_SEM_81Q_ant_num_pm{3}=std(ant_num_5um{4}./25,'omitnan')./(sqrt(length(ant_num_5um{4})));
lb_SEM_81Q_ant_num_pm{4}=std(ant_num_10um{4}./50,'omitnan')./(sqrt(length(ant_num_10um{4})));
lb_SEM_81Q_ant_num_pm{5}=std(ant_num_15um{4}./75,'omitnan')./(sqrt(length(ant_num_15um{4})));

% all_SEM_81Q_ant_num{1}=std(ant_num_pt25um{4},'omitnan')./(sqrt(length(ant_num_pt25um{4})));
% all_SEM_81Q_ant_num{2}=std(ant_num_pt5um{4},'omitnan')./(sqrt(length(ant_num_pt5um{4})));
% all_SEM_81Q_ant_num{3}=std(ant_num_1um{4},'omitnan')./(sqrt(length(ant_num_1um{4})));
% all_SEM_81Q_ant_num{4}=std(ant_num_1pt5um{4},'omitnan')./(sqrt(length(ant_num_1pt5um{4})));
% all_SEM_81Q_ant_num{5}=std(ant_num_2um{4},'omitnan')./(sqrt(length(ant_num_2um{4})));
% all_SEM_81Q_ant_num{6}=std(ant_num_5um{4},'omitnan')./(sqrt(length(ant_num_5um{4})));
% all_SEM_81Q_ant_num{7}=std(ant_num_10um{4},'omitnan')./(sqrt(length(ant_num_10um{4})));
% all_SEM_81Q_ant_num{8}=std(ant_num_15um{4},'omitnan')./(sqrt(length(ant_num_15um{4})));

% %Retrograde Fraction Means
% all_mean_30Q_ret_frac=[mean(frac_ret_1um{1},'omitnan'), mean(frac_ret_2um{1},'omitnan'),mean(frac_ret_5um{1},'omitnan'),mean(frac_ret_10um{1},'omitnan'),mean(frac_ret_15um{1},'omitnan')];
% all_mean_45Q_ret_frac=[mean(frac_ret_1um{2},'omitnan'), mean(frac_ret_2um{2},'omitnan'),mean(frac_ret_5um{2},'omitnan'),mean(frac_ret_10um{2},'omitnan'),mean(frac_ret_15um{2},'omitnan')];
% all_mean_65Q_ret_frac=[mean(frac_ret_1um{3},'omitnan'), mean(frac_ret_2um{3},'omitnan'),mean(frac_ret_5um{3},'omitnan'),mean(frac_ret_10um{3},'omitnan'),mean(frac_ret_15um{3},'omitnan')];
% all_mean_81Q_ret_frac=[mean(frac_ret_1um{4},'omitnan'), mean(frac_ret_2um{4},'omitnan'),mean(frac_ret_5um{4},'omitnan'),mean(frac_ret_10um{4},'omitnan'),mean(frac_ret_15um{4},'omitnan')];
% 
% %Retrograde Fraction: calculating SEM for errorbars
% all_SEM_30Q_ret_frac{1}=std(frac_ret_1um{1})./(sqrt(length(frac_ret_1um{1})));
% all_SEM_30Q_ret_frac{2}=std(frac_ret_2um{1})./(sqrt(length(frac_ret_2um{1})));
% all_SEM_30Q_ret_frac{3}=std(frac_ret_5um{1})./(sqrt(length(frac_ret_5um{1})));
% all_SEM_30Q_ret_frac{4}=std(frac_ret_10um{1})./(sqrt(length(frac_ret_10um{1})));
% all_SEM_30Q_ret_frac{5}=std(frac_ret_15um{1})./(sqrt(length(frac_ret_15um{1})));
% 
% all_SEM_45Q_ret_frac{1}=std(frac_ret_1um{2})./(sqrt(length(frac_ret_1um{2})));
% all_SEM_45Q_ret_frac{2}=std(frac_ret_2um{2})./(sqrt(length(frac_ret_2um{2})));
% all_SEM_45Q_ret_frac{3}=std(frac_ret_5um{2})./(sqrt(length(frac_ret_5um{2})));
% all_SEM_45Q_ret_frac{4}=std(frac_ret_10um{2})./(sqrt(length(frac_ret_10um{2})));
% all_SEM_45Q_ret_frac{5}=std(frac_ret_15um{2})./(sqrt(length(frac_ret_15um{2})));
% 
% all_SEM_65Q_ret_frac{1}=std(frac_ret_1um{3})./(sqrt(length(frac_ret_1um{3})));
% all_SEM_65Q_ret_frac{2}=std(frac_ret_2um{3})./(sqrt(length(frac_ret_2um{3})));
% all_SEM_65Q_ret_frac{3}=std(frac_ret_5um{3})./(sqrt(length(frac_ret_5um{3})));
% all_SEM_65Q_ret_frac{4}=std(frac_ret_10um{3})./(sqrt(length(frac_ret_10um{3})));
% all_SEM_65Q_ret_frac{5}=std(frac_ret_15um{3})./(sqrt(length(frac_ret_15um{3})));
% 
% all_SEM_81Q_ret_frac{1}=std(frac_ret_1um{4})./(sqrt(length(frac_ret_1um{4})));
% all_SEM_81Q_ret_frac{2}=std(frac_ret_2um{4})./(sqrt(length(frac_ret_2um{4})));
% all_SEM_81Q_ret_frac{3}=std(frac_ret_5um{4})./(sqrt(length(frac_ret_5um{4})));
% all_SEM_81Q_ret_frac{4}=std(frac_ret_10um{4})./(sqrt(length(frac_ret_10um{4})));
% all_SEM_81Q_ret_frac{5}=std(frac_ret_15um{4})./(sqrt(length(frac_ret_15um{4})));

%Retrograde Mean Number of Cargoes
lb_mean_30Q_ret_num=[mean(ret_num_1um{1},'omitnan'), mean(ret_num_2um{1},'omitnan'),mean(ret_num_5um{1},'omitnan'),mean(ret_num_10um{1},'omitnan'),mean(ret_num_15um{1},'omitnan')];
lb_mean_45Q_ret_num=[mean(ret_num_1um{2},'omitnan'), mean(ret_num_2um{2},'omitnan'),mean(ret_num_5um{2},'omitnan'),mean(ret_num_10um{2},'omitnan'),mean(ret_num_15um{2},'omitnan')];
lb_mean_65Q_ret_num=[mean(ret_num_1um{3},'omitnan'), mean(ret_num_2um{3},'omitnan'),mean(ret_num_5um{3},'omitnan'),mean(ret_num_10um{3},'omitnan'),mean(ret_num_15um{3},'omitnan')];
lb_mean_81Q_ret_num=[mean(ret_num_1um{4},'omitnan'), mean(ret_num_2um{4},'omitnan'),mean(ret_num_5um{4},'omitnan'),mean(ret_num_10um{4},'omitnan'),mean(ret_num_15um{4},'omitnan')];

lb_mean_30Q_ret_num_pm=[mean(ret_num_1um{1}./5,'omitnan'), mean(ret_num_2um{1}./10,'omitnan'),mean(ret_num_5um{1}./25,'omitnan'),mean(ret_num_10um{1}./50,'omitnan'),mean(ret_num_15um{1}./75,'omitnan')];
lb_mean_45Q_ret_num_pm=[mean(ret_num_1um{2}./5,'omitnan'), mean(ret_num_2um{2}./10,'omitnan'),mean(ret_num_5um{2}./25,'omitnan'),mean(ret_num_10um{2}./50,'omitnan'),mean(ret_num_15um{2}./75,'omitnan')];
lb_mean_65Q_ret_num_pm=[mean(ret_num_1um{3}./5,'omitnan'), mean(ret_num_2um{3}./10,'omitnan'),mean(ret_num_5um{3}./25,'omitnan'),mean(ret_num_10um{3}./50,'omitnan'),mean(ret_num_15um{3}./75,'omitnan')];
lb_mean_81Q_ret_num_pm=[mean(ret_num_1um{4}./5,'omitnan'), mean(ret_num_2um{4}./10,'omitnan'),mean(ret_num_5um{4}./25,'omitnan'),mean(ret_num_10um{4}./50,'omitnan'),mean(ret_num_15um{4}./75,'omitnan')];


% all_mean_30Q_ret_num=[mean(ret_num_pt25um{1},'omitnan'),mean(ret_num_pt5um{1},'omitnan'),mean(ret_num_1um{1},'omitnan'),mean(ret_num_1pt5um{1},'omitnan'),mean(ret_num_2um{1},'omitnan'),mean(ret_num_5um{1},'omitnan'),mean(ret_num_10um{1},'omitnan'),mean(ret_num_15um{1},'omitnan')];
% all_mean_45Q_ret_num=[mean(ret_num_pt25um{2},'omitnan'),mean(ret_num_pt5um{2},'omitnan'),mean(ret_num_1um{2},'omitnan'),mean(ret_num_1pt5um{2},'omitnan'),mean(ret_num_2um{2},'omitnan'),mean(ret_num_5um{2},'omitnan'),mean(ret_num_10um{2},'omitnan'),mean(ret_num_15um{2},'omitnan')];
% all_mean_65Q_ret_num=[mean(ret_num_pt25um{3},'omitnan'),mean(ret_num_pt5um{3},'omitnan'),mean(ret_num_1um{3},'omitnan'),mean(ret_num_1pt5um{3},'omitnan'),mean(ret_num_2um{3},'omitnan'),mean(ret_num_5um{3},'omitnan'),mean(ret_num_10um{3},'omitnan'),mean(ret_num_15um{3},'omitnan')];
% all_mean_81Q_ret_num=[mean(ret_num_pt25um{4},'omitnan'),mean(ret_num_pt5um{4},'omitnan'),mean(ret_num_1um{4},'omitnan'),mean(ret_num_1pt5um{4},'omitnan'),mean(ret_num_2um{4},'omitnan'),mean(ret_num_5um{4},'omitnan'),mean(ret_num_10um{4},'omitnan'),mean(ret_num_15um{4},'omitnan')];

%Retrograde number: calculating SEM for errorbars
lb_SEM_30Q_ret_num{1}=std(ret_num_1um{1},'omitnan')./(sqrt(length(ret_num_1um{1})));
lb_SEM_30Q_ret_num{2}=std(ret_num_2um{1},'omitnan')./(sqrt(length(ret_num_2um{1})));
lb_SEM_30Q_ret_num{3}=std(ret_num_5um{1},'omitnan')./(sqrt(length(ret_num_5um{1})));
lb_SEM_30Q_ret_num{4}=std(ret_num_10um{1},'omitnan')./(sqrt(length(ret_num_10um{1})));
lb_SEM_30Q_ret_num{5}=std(ret_num_15um{1},'omitnan')./(sqrt(length(ret_num_15um{1})));

lb_SEM_30Q_ret_num_pm{1}=std(ret_num_1um{1}./5,'omitnan')./(sqrt(length(ret_num_1um{1})));
lb_SEM_30Q_ret_num_pm{2}=std(ret_num_2um{1}./10,'omitnan')./(sqrt(length(ret_num_2um{1})));
lb_SEM_30Q_ret_num_pm{3}=std(ret_num_5um{1}./25,'omitnan')./(sqrt(length(ret_num_5um{1})));
lb_SEM_30Q_ret_num_pm{4}=std(ret_num_10um{1}./50,'omitnan')./(sqrt(length(ret_num_10um{1})));
lb_SEM_30Q_ret_num_pm{5}=std(ret_num_15um{1}./75,'omitnan')./(sqrt(length(ret_num_15um{1})));

% all_SEM_30Q_ret_num{1}=std(ret_num_pt25um{1},'omitnan')./(sqrt(length(ret_num_pt25um{1})));
% all_SEM_30Q_ret_num{2}=std(ret_num_pt5um{1},'omitnan')./(sqrt(length(ret_num_pt5um{1})));
% all_SEM_30Q_ret_num{3}=std(ret_num_1um{1},'omitnan')./(sqrt(length(ret_num_1um{1})));
% all_SEM_30Q_ret_num{4}=std(ret_num_1pt5um{1},'omitnan')./(sqrt(length(ret_num_1pt5um{1})));
% all_SEM_30Q_ret_num{5}=std(ret_num_2um{1},'omitnan')./(sqrt(length(ret_num_2um{1})));
% all_SEM_30Q_ret_num{6}=std(ret_num_5um{1},'omitnan')./(sqrt(length(ret_num_5um{1})));
% all_SEM_30Q_ret_num{7}=std(ret_num_10um{1},'omitnan')./(sqrt(length(ret_num_10um{1})));
% all_SEM_30Q_ret_num{8}=std(ret_num_15um{1},'omitnan')./(sqrt(length(ret_num_15um{1})));

lb_SEM_45Q_ret_num{1}=std(ret_num_1um{2},'omitnan')./(sqrt(length(ret_num_1um{2})));
lb_SEM_45Q_ret_num{2}=std(ret_num_2um{2},'omitnan')./(sqrt(length(ret_num_2um{2})));
lb_SEM_45Q_ret_num{3}=std(ret_num_5um{2},'omitnan')./(sqrt(length(ret_num_5um{2})));
lb_SEM_45Q_ret_num{4}=std(ret_num_10um{2},'omitnan')./(sqrt(length(ret_num_10um{2})));
lb_SEM_45Q_ret_num{5}=std(ret_num_15um{2},'omitnan')./(sqrt(length(ret_num_15um{2})));

lb_SEM_45Q_ret_num_pm{1}=std(ret_num_1um{2}./5,'omitnan')./(sqrt(length(ret_num_1um{2})));
lb_SEM_45Q_ret_num_pm{2}=std(ret_num_2um{2}./10,'omitnan')./(sqrt(length(ret_num_2um{2})));
lb_SEM_45Q_ret_num_pm{3}=std(ret_num_5um{2}./25,'omitnan')./(sqrt(length(ret_num_5um{2})));
lb_SEM_45Q_ret_num_pm{4}=std(ret_num_10um{2}./50,'omitnan')./(sqrt(length(ret_num_10um{2})));
lb_SEM_45Q_ret_num_pm{5}=std(ret_num_15um{2}./75,'omitnan')./(sqrt(length(ret_num_15um{2})));

% all_SEM_45Q_ret_num{1}=std(ret_num_pt25um{2},'omitnan')./(sqrt(length(ret_num_pt25um{2})));
% all_SEM_45Q_ret_num{2}=std(ret_num_pt5um{2},'omitnan')./(sqrt(length(ret_num_pt5um{2})));
% all_SEM_45Q_ret_num{3}=std(ret_num_1um{2},'omitnan')./(sqrt(length(ret_num_1um{2})));
% all_SEM_45Q_ret_num{4}=std(ret_num_1pt5um{2},'omitnan')./(sqrt(length(ret_num_1pt5um{2})));
% all_SEM_45Q_ret_num{5}=std(ret_num_2um{2},'omitnan')./(sqrt(length(ret_num_2um{2})));
% all_SEM_45Q_ret_num{6}=std(ret_num_5um{2},'omitnan')./(sqrt(length(ret_num_5um{2})));
% all_SEM_45Q_ret_num{7}=std(ret_num_10um{2},'omitnan')./(sqrt(length(ret_num_10um{2})));
% all_SEM_45Q_ret_num{8}=std(ret_num_15um{2},'omitnan')./(sqrt(length(ret_num_15um{2})));

lb_SEM_65Q_ret_num{1}=std(ret_num_1um{3},'omitnan')./(sqrt(length(ret_num_1um{3})));
lb_SEM_65Q_ret_num{2}=std(ret_num_2um{3},'omitnan')./(sqrt(length(ret_num_2um{3})));
lb_SEM_65Q_ret_num{3}=std(ret_num_5um{3},'omitnan')./(sqrt(length(ret_num_5um{3})));
lb_SEM_65Q_ret_num{4}=std(ret_num_10um{3},'omitnan')./(sqrt(length(ret_num_10um{3})));
lb_SEM_65Q_ret_num{5}=std(ret_num_15um{3},'omitnan')./(sqrt(length(ret_num_15um{3})));

lb_SEM_65Q_ret_num_pm{1}=std(ret_num_1um{3}./5,'omitnan')./(sqrt(length(ret_num_1um{3})));
lb_SEM_65Q_ret_num_pm{2}=std(ret_num_2um{3}./10,'omitnan')./(sqrt(length(ret_num_2um{3})));
lb_SEM_65Q_ret_num_pm{3}=std(ret_num_5um{3}./25,'omitnan')./(sqrt(length(ret_num_5um{3})));
lb_SEM_65Q_ret_num_pm{4}=std(ret_num_10um{3}./50,'omitnan')./(sqrt(length(ret_num_10um{3})));
lb_SEM_65Q_ret_num_pm{5}=std(ret_num_15um{3}./75,'omitnan')./(sqrt(length(ret_num_15um{3})));

% all_SEM_65Q_ret_num{1}=std(ret_num_pt25um{3},'omitnan')./(sqrt(length(ret_num_pt25um{3})));
% all_SEM_65Q_ret_num{2}=std(ret_num_pt5um{3},'omitnan')./(sqrt(length(ret_num_pt5um{3})));
% all_SEM_65Q_ret_num{3}=std(ret_num_1um{3},'omitnan')./(sqrt(length(ret_num_1um{3})));
% all_SEM_65Q_ret_num{4}=std(ret_num_1pt5um{3},'omitnan')./(sqrt(length(ret_num_1pt5um{3})));
% all_SEM_65Q_ret_num{5}=std(ret_num_2um{3},'omitnan')./(sqrt(length(ret_num_2um{3})));
% all_SEM_65Q_ret_num{6}=std(ret_num_5um{3},'omitnan')./(sqrt(length(ret_num_5um{3})));
% all_SEM_65Q_ret_num{7}=std(ret_num_10um{3},'omitnan')./(sqrt(length(ret_num_10um{3})));
% all_SEM_65Q_ret_num{8}=std(ret_num_15um{3},'omitnan')./(sqrt(length(ret_num_15um{3})));

lb_SEM_81Q_ret_num{1}=std(ret_num_1um{4},'omitnan')./(sqrt(length(ret_num_1um{4})));
lb_SEM_81Q_ret_num{2}=std(ret_num_2um{4},'omitnan')./(sqrt(length(ret_num_2um{4})));
lb_SEM_81Q_ret_num{3}=std(ret_num_5um{4},'omitnan')./(sqrt(length(ret_num_5um{4})));
lb_SEM_81Q_ret_num{4}=std(ret_num_10um{4},'omitnan')./(sqrt(length(ret_num_10um{4})));
lb_SEM_81Q_ret_num{5}=std(ret_num_15um{4},'omitnan')./(sqrt(length(ret_num_15um{4})));

lb_SEM_81Q_ret_num_pm{1}=std(ret_num_1um{4}./5,'omitnan')./(sqrt(length(ret_num_1um{4})));
lb_SEM_81Q_ret_num_pm{2}=std(ret_num_2um{4}./10,'omitnan')./(sqrt(length(ret_num_2um{4})));
lb_SEM_81Q_ret_num_pm{3}=std(ret_num_5um{4}./25,'omitnan')./(sqrt(length(ret_num_5um{4})));
lb_SEM_81Q_ret_num_pm{4}=std(ret_num_10um{4}./50,'omitnan')./(sqrt(length(ret_num_10um{4})));
lb_SEM_81Q_ret_num_pm{5}=std(ret_num_15um{4}./75,'omitnan')./(sqrt(length(ret_num_15um{4})));

% all_SEM_81Q_ret_num{1}=std(ret_num_pt25um{4},'omitnan')./(sqrt(length(ret_num_pt25um{4})));
% all_SEM_81Q_ret_num{2}=std(ret_num_pt5um{4},'omitnan')./(sqrt(length(ret_num_pt5um{4})));
% all_SEM_81Q_ret_num{3}=std(ret_num_1um{4},'omitnan')./(sqrt(length(ret_num_1um{4})));
% all_SEM_81Q_ret_num{4}=std(ret_num_1pt5um{4},'omitnan')./(sqrt(length(ret_num_1pt5um{4})));
% all_SEM_81Q_ret_num{5}=std(ret_num_2um{4},'omitnan')./(sqrt(length(ret_num_2um{4})));
% all_SEM_81Q_ret_num{6}=std(ret_num_5um{4},'omitnan')./(sqrt(length(ret_num_5um{4})));
% all_SEM_81Q_ret_num{7}=std(ret_num_10um{4},'omitnan')./(sqrt(length(ret_num_10um{4})));
% all_SEM_81Q_ret_num{8}=std(ret_num_15um{4},'omitnan')./(sqrt(length(ret_num_15um{4})));

%No net flux Fraction Means
all_mean_30Q_no_net_frac=[mean(frac_no_net_1um{1},'omitnan'), mean(frac_no_net_2um{1},'omitnan'),mean(frac_no_net_5um{1},'omitnan'),mean(frac_no_net_10um{1},'omitnan'),mean(frac_no_net_15um{1},'omitnan')];
all_mean_45Q_no_net_frac=[mean(frac_no_net_1um{2},'omitnan'), mean(frac_no_net_2um{2},'omitnan'),mean(frac_no_net_5um{2},'omitnan'),mean(frac_no_net_10um{2},'omitnan'),mean(frac_no_net_15um{2},'omitnan')];
all_mean_65Q_no_net_frac=[mean(frac_no_net_1um{3},'omitnan'), mean(frac_no_net_2um{3},'omitnan'),mean(frac_no_net_5um{3},'omitnan'),mean(frac_no_net_10um{3},'omitnan'),mean(frac_no_net_15um{3},'omitnan')];
all_mean_81Q_no_net_frac=[mean(frac_no_net_1um{4},'omitnan'), mean(frac_no_net_2um{4},'omitnan'),mean(frac_no_net_5um{4},'omitnan'),mean(frac_no_net_10um{4},'omitnan'),mean(frac_no_net_15um{4},'omitnan')];

lb_mean_30Q_no_net_frac_pm=[mean(frac_no_net_1um{1}./5,'omitnan'), mean(frac_no_net_2um{1}./10,'omitnan'),mean(frac_no_net_5um{1}./25,'omitnan'),mean(frac_no_net_10um{1}./50,'omitnan'),mean(frac_no_net_15um{1}./75,'omitnan')];
lb_mean_45Q_no_net_frac_pm=[mean(frac_no_net_1um{2}./5,'omitnan'), mean(frac_no_net_2um{2}./10,'omitnan'),mean(frac_no_net_5um{2}./25,'omitnan'),mean(frac_no_net_10um{2}./50,'omitnan'),mean(frac_no_net_15um{2}./75,'omitnan')];
lb_mean_65Q_no_net_frac_pm=[mean(frac_no_net_1um{3}./5,'omitnan'), mean(frac_no_net_2um{3}./10,'omitnan'),mean(frac_no_net_5um{3}./25,'omitnan'),mean(frac_no_net_10um{3}./50,'omitnan'),mean(frac_no_net_15um{3}./75,'omitnan')];
lb_mean_81Q_no_net_frac_pm=[mean(frac_no_net_1um{4}./5,'omitnan'), mean(frac_no_net_2um{4}./10,'omitnan'),mean(frac_no_net_5um{4}./25,'omitnan'),mean(frac_no_net_10um{4}./50,'omitnan'),mean(frac_no_net_15um{4}./75,'omitnan')];

%No net flux Fraction: calculating SEM for errorbars
all_SEM_30Q_no_net_frac{1}=std(frac_no_net_1um{1},'omitnan')./(sqrt(length(frac_no_net_1um{1})));
all_SEM_30Q_no_net_frac{2}=std(frac_no_net_2um{1},'omitnan')./(sqrt(length(frac_no_net_2um{1})));
all_SEM_30Q_no_net_frac{3}=std(frac_no_net_5um{1},'omitnan')./(sqrt(length(frac_no_net_5um{1})));
all_SEM_30Q_no_net_frac{4}=std(frac_no_net_10um{1},'omitnan')./(sqrt(length(frac_no_net_10um{1})));
all_SEM_30Q_no_net_frac{5}=std(frac_no_net_15um{1},'omitnan')./(sqrt(length(frac_no_net_15um{1})));

all_SEM_30Q_no_net_frac_pm{1}=std(frac_no_net_1um{1}./5,'omitnan')./(sqrt(length(frac_no_net_1um{1})));
all_SEM_30Q_no_net_frac_pm{2}=std(frac_no_net_2um{1}./10,'omitnan')./(sqrt(length(frac_no_net_2um{1})));
all_SEM_30Q_no_net_frac_pm{3}=std(frac_no_net_5um{1}./25,'omitnan')./(sqrt(length(frac_no_net_5um{1})));
all_SEM_30Q_no_net_frac_pm{4}=std(frac_no_net_10um{1}./50,'omitnan')./(sqrt(length(frac_no_net_10um{1})));
all_SEM_30Q_no_net_frac_pm{5}=std(frac_no_net_15um{1}./75,'omitnan')./(sqrt(length(frac_no_net_15um{1})));

all_SEM_45Q_no_net_frac{1}=std(frac_no_net_1um{2},'omitnan')./(sqrt(length(frac_no_net_1um{2})));
all_SEM_45Q_no_net_frac{2}=std(frac_no_net_2um{2},'omitnan')./(sqrt(length(frac_no_net_2um{2})));
all_SEM_45Q_no_net_frac{3}=std(frac_no_net_5um{2},'omitnan')./(sqrt(length(frac_no_net_5um{2})));
all_SEM_45Q_no_net_frac{4}=std(frac_no_net_10um{2},'omitnan')./(sqrt(length(frac_no_net_10um{2})));
all_SEM_45Q_no_net_frac{5}=std(frac_no_net_15um{2},'omitnan')./(sqrt(length(frac_no_net_15um{2})));

all_SEM_45Q_no_net_frac_pm{1}=std(frac_no_net_1um{2},'omitnan')./(sqrt(length(frac_no_net_1um{2})));
all_SEM_45Q_no_net_frac_pm{2}=std(frac_no_net_2um{2},'omitnan')./(sqrt(length(frac_no_net_2um{2})));
all_SEM_45Q_no_net_frac_pm{3}=std(frac_no_net_5um{2},'omitnan')./(sqrt(length(frac_no_net_5um{2})));
all_SEM_45Q_no_net_frac_pm{4}=std(frac_no_net_10um{2},'omitnan')./(sqrt(length(frac_no_net_10um{2})));
all_SEM_45Q_no_net_frac_pm{5}=std(frac_no_net_15um{2},'omitnan')./(sqrt(length(frac_no_net_15um{2})));

all_SEM_65Q_no_net_frac{1}=std(frac_no_net_1um{3},'omitnan')./(sqrt(length(frac_no_net_1um{3})));
all_SEM_65Q_no_net_frac{2}=std(frac_no_net_2um{3},'omitnan')./(sqrt(length(frac_no_net_2um{3})));
all_SEM_65Q_no_net_frac{3}=std(frac_no_net_5um{3},'omitnan')./(sqrt(length(frac_no_net_5um{3})));
all_SEM_65Q_no_net_frac{4}=std(frac_no_net_10um{3},'omitnan')./(sqrt(length(frac_no_net_10um{3})));
all_SEM_65Q_no_net_frac{5}=std(frac_no_net_15um{3},'omitnan')./(sqrt(length(frac_no_net_15um{3})));

all_SEM_65Q_no_net_frac_pm{1}=std(frac_no_net_1um{3}./5,'omitnan')./(sqrt(length(frac_no_net_1um{3})));
all_SEM_65Q_no_net_frac_pm{2}=std(frac_no_net_2um{3}./10,'omitnan')./(sqrt(length(frac_no_net_2um{3})));
all_SEM_65Q_no_net_frac_pm{3}=std(frac_no_net_5um{3}./25,'omitnan')./(sqrt(length(frac_no_net_5um{3})));
all_SEM_65Q_no_net_frac_pm{4}=std(frac_no_net_10um{3}./50,'omitnan')./(sqrt(length(frac_no_net_10um{3})));
all_SEM_65Q_no_net_frac_pm{5}=std(frac_no_net_15um{3}./75,'omitnan')./(sqrt(length(frac_no_net_15um{3})));

all_SEM_81Q_no_net_frac{1}=std(frac_no_net_1um{4},'omitnan')./(sqrt(length(frac_no_net_1um{4})));
all_SEM_81Q_no_net_frac{2}=std(frac_no_net_2um{4},'omitnan')./(sqrt(length(frac_no_net_2um{4})));
all_SEM_81Q_no_net_frac{3}=std(frac_no_net_5um{4},'omitnan')./(sqrt(length(frac_no_net_5um{4})));
all_SEM_81Q_no_net_frac{4}=std(frac_no_net_10um{4},'omitnan')./(sqrt(length(frac_no_net_10um{4})));
all_SEM_81Q_no_net_frac{5}=std(frac_no_net_15um{4},'omitnan')./(sqrt(length(frac_no_net_15um{4})));

all_SEM_81Q_no_net_frac_pm{1}=std(frac_no_net_1um{4}./5,'omitnan')./(sqrt(length(frac_no_net_1um{4})));
all_SEM_81Q_no_net_frac_pm{2}=std(frac_no_net_2um{4}./10,'omitnan')./(sqrt(length(frac_no_net_2um{4})));
all_SEM_81Q_no_net_frac_pm{3}=std(frac_no_net_5um{4}./25,'omitnan')./(sqrt(length(frac_no_net_5um{4})));
all_SEM_81Q_no_net_frac_pm{4}=std(frac_no_net_10um{4}./50,'omitnan')./(sqrt(length(frac_no_net_10um{4})));
all_SEM_81Q_no_net_frac_pm{5}=std(frac_no_net_15um{4}./75,'omitnan')./(sqrt(length(frac_no_net_15um{4})));

%No net flux Mean Number of Cargoes
lb_mean_30Q_no_net_num=[mean(no_net_num_1um{1},'omitnan'), mean(no_net_num_2um{1},'omitnan'),mean(no_net_num_5um{1},'omitnan'),mean(no_net_num_10um{1},'omitnan'),mean(no_net_num_15um{1},'omitnan')];
lb_mean_45Q_no_net_num=[mean(no_net_num_1um{2},'omitnan'), mean(no_net_num_2um{2},'omitnan'),mean(no_net_num_5um{2},'omitnan'),mean(no_net_num_10um{2},'omitnan'),mean(no_net_num_15um{2},'omitnan')];
lb_mean_65Q_no_net_num=[mean(no_net_num_1um{3},'omitnan'), mean(no_net_num_2um{3},'omitnan'),mean(no_net_num_5um{3},'omitnan'),mean(no_net_num_10um{3},'omitnan'),mean(no_net_num_15um{3},'omitnan')];
lb_mean_81Q_no_net_num=[mean(no_net_num_1um{4},'omitnan'), mean(no_net_num_2um{4},'omitnan'),mean(no_net_num_5um{4},'omitnan'),mean(no_net_num_10um{4},'omitnan'),mean(no_net_num_15um{4},'omitnan')];

lb_mean_30Q_no_net_num_pm=[mean(no_net_num_1um{1}./5,'omitnan'), mean(no_net_num_2um{1}./10,'omitnan'),mean(no_net_num_5um{1}./25,'omitnan'),mean(no_net_num_10um{1}./50,'omitnan'),mean(no_net_num_15um{1}./75,'omitnan')];
lb_mean_45Q_no_net_num_pm=[mean(no_net_num_1um{2}./5,'omitnan'), mean(no_net_num_2um{2}./10,'omitnan'),mean(no_net_num_5um{2}./25,'omitnan'),mean(no_net_num_10um{2}./50,'omitnan'),mean(no_net_num_15um{2}./75,'omitnan')];
lb_mean_65Q_no_net_num_pm=[mean(no_net_num_1um{3}./5,'omitnan'), mean(no_net_num_2um{3}./10,'omitnan'),mean(no_net_num_5um{3}./25,'omitnan'),mean(no_net_num_10um{3}./50,'omitnan'),mean(no_net_num_15um{3}./75,'omitnan')];
lb_mean_81Q_no_net_num_pm=[mean(no_net_num_1um{4}./5,'omitnan'), mean(no_net_num_2um{4}./10,'omitnan'),mean(no_net_num_5um{4}./25,'omitnan'),mean(no_net_num_10um{4}./50,'omitnan'),mean(no_net_num_15um{4}./75,'omitnan')];

% all_mean_30Q_no_net_num=[mean(no_net_num_pt25um{1},'omitnan'),mean(no_net_num_pt5um{1},'omitnan'),mean(no_net_num_1um{1},'omitnan'),mean(no_net_num_1pt5um{1},'omitnan'),mean(no_net_num_2um{1},'omitnan'),mean(no_net_num_5um{1},'omitnan'),mean(no_net_num_10um{1},'omitnan'),mean(no_net_num_15um{1},'omitnan')];
% all_mean_45Q_no_net_num=[mean(no_net_num_pt25um{2},'omitnan'),mean(no_net_num_pt5um{2},'omitnan'),mean(no_net_num_1um{2},'omitnan'),mean(no_net_num_1pt5um{2},'omitnan'),mean(no_net_num_2um{2},'omitnan'),mean(no_net_num_5um{2},'omitnan'),mean(no_net_num_10um{2},'omitnan'),mean(no_net_num_15um{2},'omitnan')];
% all_mean_65Q_no_net_num=[mean(no_net_num_pt25um{3},'omitnan'),mean(no_net_num_pt5um{3},'omitnan'),mean(no_net_num_1um{3},'omitnan'),mean(no_net_num_1pt5um{3},'omitnan'),mean(no_net_num_2um{3},'omitnan'),mean(no_net_num_5um{3},'omitnan'),mean(no_net_num_10um{3},'omitnan'),mean(no_net_num_15um{3},'omitnan')];
% all_mean_81Q_no_net_num=[mean(no_net_num_pt25um{4},'omitnan'),mean(no_net_num_pt5um{4},'omitnan'),mean(no_net_num_1um{4},'omitnan'),mean(no_net_num_1pt5um{4},'omitnan'),mean(no_net_num_2um{4},'omitnan'),mean(no_net_num_5um{4},'omitnan'),mean(no_net_num_10um{4},'omitnan'),mean(no_net_num_15um{4},'omitnan')];
%No net number: calculating SEM for errorbars
lb_SEM_30Q_no_net_num{1}=std(no_net_num_1um{1},'omitnan')./(sqrt(length(no_net_num_1um{1})));
lb_SEM_30Q_no_net_num{2}=std(no_net_num_2um{1},'omitnan')./(sqrt(length(no_net_num_2um{1})));
lb_SEM_30Q_no_net_num{3}=std(no_net_num_5um{1},'omitnan')./(sqrt(length(no_net_num_5um{1})));
lb_SEM_30Q_no_net_num{4}=std(no_net_num_10um{1},'omitnan')./(sqrt(length(no_net_num_10um{1})));
lb_SEM_30Q_no_net_num{5}=std(no_net_num_15um{1},'omitnan')./(sqrt(length(no_net_num_15um{1})));

lb_SEM_30Q_no_net_num_pm{1}=std(no_net_num_1um{1}./5,'omitnan')./(sqrt(length(no_net_num_1um{1})));
lb_SEM_30Q_no_net_num_pm{2}=std(no_net_num_2um{1}./10,'omitnan')./(sqrt(length(no_net_num_2um{1})));
lb_SEM_30Q_no_net_num_pm{3}=std(no_net_num_5um{1}./25,'omitnan')./(sqrt(length(no_net_num_5um{1})));
lb_SEM_30Q_no_net_num_pm{4}=std(no_net_num_10um{1}./50,'omitnan')./(sqrt(length(no_net_num_10um{1})));
lb_SEM_30Q_no_net_num_pm{5}=std(no_net_num_15um{1}./75,'omitnan')./(sqrt(length(no_net_num_15um{1})));


% all_SEM_30Q_no_net_num{1}=std(no_net_num_pt25um{1},'omitnan')./(sqrt(length(no_net_num_pt25um{1})));
% all_SEM_30Q_no_net_num{2}=std(no_net_num_pt5um{1},'omitnan')./(sqrt(length(no_net_num_pt5um{1})));
% all_SEM_30Q_no_net_num{3}=std(no_net_num_1um{1},'omitnan')./(sqrt(length(no_net_num_1um{1})));
% all_SEM_30Q_no_net_num{4}=std(no_net_num_1pt5um{1},'omitnan')./(sqrt(length(no_net_num_1pt5um{1})));
% all_SEM_30Q_no_net_num{5}=std(no_net_num_2um{1},'omitnan')./(sqrt(length(no_net_num_2um{1})));
% all_SEM_30Q_no_net_num{6}=std(no_net_num_5um{1},'omitnan')./(sqrt(length(no_net_num_5um{1})));
% all_SEM_30Q_no_net_num{7}=std(no_net_num_10um{1},'omitnan')./(sqrt(length(no_net_num_10um{1})));
% all_SEM_30Q_no_net_num{8}=std(no_net_num_15um{1},'omitnan')./(sqrt(length(no_net_num_15um{1})));

lb_SEM_45Q_no_net_num{1}=std(no_net_num_1um{2},'omitnan')./(sqrt(length(no_net_num_1um{2})));
lb_SEM_45Q_no_net_num{2}=std(no_net_num_2um{2},'omitnan')./(sqrt(length(no_net_num_2um{2})));
lb_SEM_45Q_no_net_num{3}=std(no_net_num_5um{2},'omitnan')./(sqrt(length(no_net_num_5um{2})));
lb_SEM_45Q_no_net_num{4}=std(no_net_num_10um{2},'omitnan')./(sqrt(length(no_net_num_10um{2})));
lb_SEM_45Q_no_net_num{5}=std(no_net_num_15um{2},'omitnan')./(sqrt(length(no_net_num_15um{2})));

lb_SEM_45Q_no_net_num_pm{1}=std(no_net_num_1um{2}./5,'omitnan')./(sqrt(length(no_net_num_1um{2})));
lb_SEM_45Q_no_net_num_pm{2}=std(no_net_num_2um{2}./10,'omitnan')./(sqrt(length(no_net_num_2um{2})));
lb_SEM_45Q_no_net_num_pm{3}=std(no_net_num_5um{2}./25,'omitnan')./(sqrt(length(no_net_num_5um{2})));
lb_SEM_45Q_no_net_num_pm{4}=std(no_net_num_10um{2}./50,'omitnan')./(sqrt(length(no_net_num_10um{2})));
lb_SEM_45Q_no_net_num_pm{5}=std(no_net_num_15um{2}./75,'omitnan')./(sqrt(length(no_net_num_15um{2})));

% all_SEM_45Q_no_net_num{1}=std(no_net_num_pt25um{2},'omitnan')./(sqrt(length(no_net_num_pt25um{2})));
% all_SEM_45Q_no_net_num{2}=std(no_net_num_pt5um{2},'omitnan')./(sqrt(length(no_net_num_pt5um{2})));
% all_SEM_45Q_no_net_num{3}=std(no_net_num_1um{2},'omitnan')./(sqrt(length(no_net_num_1um{2})));
% all_SEM_45Q_no_net_num{4}=std(no_net_num_1pt5um{2},'omitnan')./(sqrt(length(no_net_num_1pt5um{2})));
% all_SEM_45Q_no_net_num{5}=std(no_net_num_2um{2},'omitnan')./(sqrt(length(no_net_num_2um{2})));
% all_SEM_45Q_no_net_num{6}=std(no_net_num_5um{2},'omitnan')./(sqrt(length(no_net_num_5um{2})));
% all_SEM_45Q_no_net_num{7}=std(no_net_num_10um{2},'omitnan')./(sqrt(length(no_net_num_10um{2})));
% all_SEM_45Q_no_net_num{8}=std(no_net_num_15um{2},'omitnan')./(sqrt(length(no_net_num_15um{2})));

lb_SEM_65Q_no_net_num{1}=std(no_net_num_1um{3},'omitnan')./(sqrt(length(no_net_num_1um{3})));
lb_SEM_65Q_no_net_num{2}=std(no_net_num_2um{3},'omitnan')./(sqrt(length(no_net_num_2um{3})));
lb_SEM_65Q_no_net_num{3}=std(no_net_num_5um{3},'omitnan')./(sqrt(length(no_net_num_5um{3})));
lb_SEM_65Q_no_net_num{4}=std(no_net_num_10um{3},'omitnan')./(sqrt(length(no_net_num_10um{3})));
lb_SEM_65Q_no_net_num{5}=std(no_net_num_15um{3},'omitnan')./(sqrt(length(no_net_num_15um{3})));

lb_SEM_65Q_no_net_num_pm{1}=std(no_net_num_1um{3}./5,'omitnan')./(sqrt(length(no_net_num_1um{3})));
lb_SEM_65Q_no_net_num_pm{2}=std(no_net_num_2um{3}./10,'omitnan')./(sqrt(length(no_net_num_2um{3})));
lb_SEM_65Q_no_net_num_pm{3}=std(no_net_num_5um{3}./25,'omitnan')./(sqrt(length(no_net_num_5um{3})));
lb_SEM_65Q_no_net_num_pm{4}=std(no_net_num_10um{3}./50,'omitnan')./(sqrt(length(no_net_num_10um{3})));
lb_SEM_65Q_no_net_num_pm{5}=std(no_net_num_15um{3}./75,'omitnan')./(sqrt(length(no_net_num_15um{3})));

% all_SEM_65Q_no_net_num{1}=std(no_net_num_pt25um{3},'omitnan')./(sqrt(length(no_net_num_pt25um{3})));
% all_SEM_65Q_no_net_num{2}=std(no_net_num_pt5um{3},'omitnan')./(sqrt(length(no_net_num_pt5um{3})));
% all_SEM_65Q_no_net_num{3}=std(no_net_num_1um{3},'omitnan')./(sqrt(length(no_net_num_1um{3})));
% all_SEM_65Q_no_net_num{4}=std(no_net_num_1pt5um{3},'omitnan')./(sqrt(length(no_net_num_1pt5um{3})));
% all_SEM_65Q_no_net_num{5}=std(no_net_num_2um{3},'omitnan')./(sqrt(length(no_net_num_2um{3})));
% all_SEM_65Q_no_net_num{6}=std(no_net_num_5um{3},'omitnan')./(sqrt(length(no_net_num_5um{3})));
% all_SEM_65Q_no_net_num{7}=std(no_net_num_10um{3},'omitnan')./(sqrt(length(no_net_num_10um{3})));
% all_SEM_65Q_no_net_num{8}=std(no_net_num_15um{3},'omitnan')./(sqrt(length(no_net_num_15um{3})));

lb_SEM_81Q_no_net_num{1}=std(no_net_num_1um{4},'omitnan')./(sqrt(length(no_net_num_1um{4})));
lb_SEM_81Q_no_net_num{2}=std(no_net_num_2um{4},'omitnan')./(sqrt(length(no_net_num_2um{4})));
lb_SEM_81Q_no_net_num{3}=std(no_net_num_5um{4},'omitnan')./(sqrt(length(no_net_num_5um{4})));
lb_SEM_81Q_no_net_num{4}=std(no_net_num_10um{4},'omitnan')./(sqrt(length(no_net_num_10um{4})));
lb_SEM_81Q_no_net_num{5}=std(no_net_num_15um{4},'omitnan')./(sqrt(length(no_net_num_15um{4})));

lb_SEM_81Q_no_net_num_pm{1}=std(no_net_num_1um{4}./5,'omitnan')./(sqrt(length(no_net_num_1um{4})));
lb_SEM_81Q_no_net_num_pm{2}=std(no_net_num_2um{4}./10,'omitnan')./(sqrt(length(no_net_num_2um{4})));
lb_SEM_81Q_no_net_num_pm{3}=std(no_net_num_5um{4}./25,'omitnan')./(sqrt(length(no_net_num_5um{4})));
lb_SEM_81Q_no_net_num_pm{4}=std(no_net_num_10um{4}./50,'omitnan')./(sqrt(length(no_net_num_10um{4})));
lb_SEM_81Q_no_net_num_pm{5}=std(no_net_num_15um{4}./75,'omitnan')./(sqrt(length(no_net_num_15um{4})));

% all_SEM_81Q_no_net_num{1}=std(no_net_num_pt25um{4},'omitnan')./(sqrt(length(no_net_num_pt25um{4})));
% all_SEM_81Q_no_net_num{2}=std(no_net_num_pt5um{4},'omitnan')./(sqrt(length(no_net_num_pt5um{4})));
% all_SEM_81Q_no_net_num{3}=std(no_net_num_1um{4},'omitnan')./(sqrt(length(no_net_num_1um{4})));
% all_SEM_81Q_no_net_num{4}=std(no_net_num_1pt5um{4},'omitnan')./(sqrt(length(no_net_num_1pt5um{4})));
% all_SEM_81Q_no_net_num{5}=std(no_net_num_2um{4},'omitnan')./(sqrt(length(no_net_num_2um{4})));
% all_SEM_81Q_no_net_num{6}=std(no_net_num_5um{4},'omitnan')./(sqrt(length(no_net_num_5um{4})));
% all_SEM_81Q_no_net_num{7}=std(no_net_num_10um{4},'omitnan')./(sqrt(length(no_net_num_10um{4})));
% all_SEM_81Q_no_net_num{8}=std(no_net_num_15um{4},'omitnan')./(sqrt(length(no_net_num_15um{4})));

large_box_sizes=[1,2,5,10,15];
all_box_sizes=[0.25, 0.5, 1, 1.5, 2, 5, 10, 15];

%Plotting the data
figure('Name','tot_enters','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_tot_enters,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_tot_enters,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_tot_enters,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_tot_enters,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.1,lb_mean_30Q_tot_enters,cell2mat(lb_SEM_30Q_tot_enters),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.05,lb_mean_45Q_tot_enters,cell2mat(lb_SEM_45Q_tot_enters),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.05,lb_mean_65Q_tot_enters,cell2mat(lb_SEM_65Q_tot_enters),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_81Q_tot_enters,cell2mat(lb_SEM_81Q_tot_enters),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Total Cargoes');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

figure('Name','tot_enters_pm','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_tot_enters,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_tot_enters,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_tot_enters,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_tot_enters,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.1,lb_mean_30Q_tot_enters_pm,cell2mat(lb_SEM_30Q_tot_enters_pm),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.05,lb_mean_45Q_tot_enters_pm,cell2mat(lb_SEM_45Q_tot_enters_pm),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.05,lb_mean_65Q_tot_enters_pm,cell2mat(lb_SEM_65Q_tot_enters_pm),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_81Q_tot_enters_pm,cell2mat(lb_SEM_81Q_tot_enters_pm),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Cargoes/\mum');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

% figure('Name','tot_enters_all','NumberTitle','off'), hold on,
% % plot(all_box_sizes-0.1,all_mean_30Q_tot_enters,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% % plot(all_box_sizes-0.05,all_mean_45Q_tot_enters,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% % plot(all_box_sizes+0.05,all_mean_65Q_tot_enters,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% % plot(all_box_sizes+0.1,all_mean_81Q_tot_enters,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
% errorbar(all_box_sizes-0.1,all_mean_30Q_tot_enters,cell2mat(all_SEM_30Q_tot_enters),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes-0.05,all_mean_45Q_tot_enters,cell2mat(all_SEM_45Q_tot_enters),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.05,all_mean_65Q_tot_enters,cell2mat(all_SEM_65Q_tot_enters),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.1,all_mean_81Q_tot_enters,cell2mat(all_SEM_81Q_tot_enters),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% xticks([1 2 5 10 15]);
% xlim([0 15.5]);
% xlabel('Width (\mum)');
% ylabel('Total Cargoes');
% publication_fig(0,0,1);
% box("off");
% % pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

figure('Name','ant_num','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_ant_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_ant_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_ant_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_ant_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_ant_num,cell2mat(lb_SEM_30Q_ant_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_ant_num,cell2mat(lb_SEM_45Q_ant_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_ant_num,cell2mat(lb_SEM_65Q_ant_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_ant_num,cell2mat(lb_SEM_81Q_ant_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Anterograde Cargoes');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

figure('Name','ant_num_pm','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_ant_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_ant_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_ant_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_ant_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_ant_num_pm,cell2mat(lb_SEM_30Q_ant_num_pm),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_ant_num_pm,cell2mat(lb_SEM_45Q_ant_num_pm),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_ant_num_pm,cell2mat(lb_SEM_65Q_ant_num_pm),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_ant_num_pm,cell2mat(lb_SEM_81Q_ant_num_pm),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Anterograde Cargoes/\mum');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

% figure('Name','ant_num_all','NumberTitle','off'), hold on,
% % plot(all_box_sizes-0.1,all_mean_30Q_ant_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% % plot(all_box_sizes-0.05,all_mean_45Q_ant_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% % plot(all_box_sizes+0.05,all_mean_65Q_ant_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% % plot(all_box_sizes+0.1,all_mean_81Q_ant_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
% errorbar(all_box_sizes-0.2,all_mean_30Q_ant_num,cell2mat(all_SEM_30Q_ant_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes-0.1,all_mean_45Q_ant_num,cell2mat(all_SEM_45Q_ant_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.1,all_mean_65Q_ant_num,cell2mat(all_SEM_65Q_ant_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.2,all_mean_81Q_ant_num,cell2mat(all_SEM_81Q_ant_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% xticks([1 2 5 10 15]);
% xlim([0 15.5]);
% xlabel('Width (\mum)');
% ylabel('Anterograde Cargoes');
% publication_fig(0,0,1);
% box("off");
% % pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

figure('Name','ret_num','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_ret_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_ret_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_ret_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_ret_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_ret_num,cell2mat(lb_SEM_30Q_ret_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_ret_num,cell2mat(lb_SEM_45Q_ret_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_ret_num,cell2mat(lb_SEM_65Q_ret_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_ret_num,cell2mat(lb_SEM_81Q_ret_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Retrograde Cargoes');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

figure('Name','ret_num_pm','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_ret_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_ret_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_ret_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_ret_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_ret_num_pm,cell2mat(lb_SEM_30Q_ret_num_pm),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_ret_num_pm,cell2mat(lb_SEM_45Q_ret_num_pm),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_ret_num_pm,cell2mat(lb_SEM_65Q_ret_num_pm),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_ret_num_pm,cell2mat(lb_SEM_81Q_ret_num_pm),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('Retrograde Cargoes/\mum');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);


% figure('Name','ret_num_all','NumberTitle','off'), hold on,
% % plot(all_box_sizes-0.1,all_mean_30Q_ret_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% % plot(all_box_sizes-0.05,all_mean_45Q_ret_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% % plot(all_box_sizes+0.05,all_mean_65Q_ret_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% % plot(all_box_sizes+0.1,all_mean_81Q_ret_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
% errorbar(all_box_sizes-0.2,all_mean_30Q_ret_num,cell2mat(all_SEM_30Q_ret_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes-0.1,all_mean_45Q_ret_num,cell2mat(all_SEM_45Q_ret_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.1,all_mean_65Q_ret_num,cell2mat(all_SEM_65Q_ret_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.2,all_mean_81Q_ret_num,cell2mat(all_SEM_81Q_ret_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% xticks([1 2 5 10 15]);
% xlim([0 15.5]);
% xlabel('Width (\mum)');
% ylabel('Retrograde Cargoes');
% publication_fig(0,0,1);
% box("off");
% % pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

figure('Name','no_net_num','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_no_net_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_no_net_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_no_net_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_no_net_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_no_net_num,cell2mat(lb_SEM_30Q_no_net_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_no_net_num,cell2mat(lb_SEM_45Q_no_net_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_no_net_num,cell2mat(lb_SEM_65Q_no_net_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_no_net_num,cell2mat(lb_SEM_81Q_no_net_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('No net flux Cargoes');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

figure('Name','no_net_num_pm','NumberTitle','off'), hold on,
% plot(large_box_sizes-0.1,lb_mean_30Q_no_net_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% plot(large_box_sizes-0.05,lb_mean_45Q_no_net_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% plot(large_box_sizes+0.05,lb_mean_65Q_no_net_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% plot(large_box_sizes+0.1,lb_mean_81Q_no_net_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
errorbar(large_box_sizes-0.2,lb_mean_30Q_no_net_num_pm,cell2mat(lb_SEM_30Q_no_net_num_pm),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes-0.1,lb_mean_45Q_no_net_num_pm,cell2mat(lb_SEM_45Q_no_net_num_pm),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.1,lb_mean_65Q_no_net_num_pm,cell2mat(lb_SEM_65Q_no_net_num_pm),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
errorbar(large_box_sizes+0.2,lb_mean_81Q_no_net_num_pm,cell2mat(lb_SEM_81Q_no_net_num_pm),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
xticks([1 2 5 10 15]);
xlim([0 15.5]);
xlabel('Width (\mum)');
ylabel('No net flux Cargoes/\mum');
publication_fig(0,0,1);
box("off");
% pbaspect([2 1 1]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
set(gca, 'FontSize',40);

% figure('Name','no_net_num_all','NumberTitle','off'), hold on,
% % plot(all_box_sizes-0.1,all_mean_30Q_no_net_num,'-','Color',colour_30Q,'Marker','.','MarkerEdgeColor',colour_30Q,'MarkerSize',20);
% % plot(all_box_sizes-0.05,all_mean_45Q_no_net_num,'-.','Color',colour_45Q,'Marker','.','MarkerEdgeColor',colour_45Q,'MarkerSize',20);
% % plot(all_box_sizes+0.05,all_mean_65Q_no_net_num,'--','Color',colour_65Q,'Marker','.','MarkerEdgeColor',colour_65Q,'MarkerSize',20);
% % plot(all_box_sizes+0.1,all_mean_81Q_no_net_num,':','Color',colour_81Q,'Marker','.','MarkerEdgeColor',colour_81Q,'MarkerSize',20);
% errorbar(all_box_sizes-0.2,all_mean_30Q_no_net_num,cell2mat(all_SEM_30Q_no_net_num),'-','Color',colour_30Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes-0.1,all_mean_45Q_no_net_num,cell2mat(all_SEM_45Q_no_net_num),'-.','Color',colour_45Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.1,all_mean_65Q_no_net_num,cell2mat(all_SEM_65Q_no_net_num),'--','Color',colour_65Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% errorbar(all_box_sizes+0.2,all_mean_81Q_no_net_num,cell2mat(all_SEM_81Q_no_net_num),':','Color',colour_81Q,'LineWidth',4,'CapSize',0,'Marker','o','MarkerSize',15);
% xticks([1 2 5 10 15]);
% xlim([0 15.5]);
% xlabel('Width (\mum)');
% ylabel('No net flux Cargoes');
% publication_fig(0,0,1);
% box("off");
% % pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Figures_motility/';  % Your destination folder
tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/';  % Your destination folder
FolderName = tempdir;  % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
 FigHandle = FigList(iFig);
 FigName  = get(FigHandle, 'Name');
 savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
 saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end

% figure('Name','ant_frac','NumberTitle','off'), hold on,
% plot(box_sizes,all_mean_30Q_ant_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,all_mean_45Q_ant_frac,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,all_mean_65Q_ant_frac,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,all_mean_81Q_ant_frac,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,all_mean_30Q_ant_frac,cell2mat(all_SEM_30Q_ant_frac),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_45Q_ant_frac,cell2mat(all_SEM_45Q_ant_frac),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_65Q_ant_frac,cell2mat(all_SEM_65Q_ant_frac),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_81Q_ant_frac,cell2mat(all_SEM_81Q_ant_frac),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Fraction of Cargoes');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

% figure('Name','ret_frac','NumberTitle','off'), hold on,
% plot(box_sizes,all_mean_30Q_ret_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,all_mean_45Q_ret_frac,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,all_mean_65Q_ret_frac,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,all_mean_81Q_ret_frac,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,all_mean_30Q_ret_frac,cell2mat(all_SEM_30Q_ret_frac),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_45Q_ret_frac,cell2mat(all_SEM_45Q_ret_frac),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_65Q_ret_frac,cell2mat(all_SEM_65Q_ret_frac),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_81Q_ret_frac,cell2mat(all_SEM_81Q_ret_frac),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Fraction of Cargoes');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);

% figure('Name','ant_ret_nums','NumberTitle','off'), hold on,
% plot(box_sizes,all_mean_30Q_ant_num,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,all_mean_45Q_ant_num,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,all_mean_65Q_ant_num,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,all_mean_81Q_ant_num,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,all_mean_30Q_ant_num,cell2mat(all_SEM_30Q_ant_num),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_45Q_ant_num,cell2mat(all_SEM_45Q_ant_num),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_65Q_ant_num,cell2mat(all_SEM_65Q_ant_num),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_81Q_ant_num,cell2mat(all_SEM_81Q_ant_num),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% plot(box_sizes,-1*all_mean_30Q_ret_num,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,-1*all_mean_45Q_ret_num,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,-1*all_mean_65Q_ret_num,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,-1*all_mean_81Q_ret_num,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,-1*all_mean_30Q_ret_num,-1*cell2mat(all_SEM_30Q_ret_num),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_45Q_ret_num,-1*cell2mat(all_SEM_45Q_ret_num),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_65Q_ret_num,-1*cell2mat(all_SEM_65Q_ret_num),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_81Q_ret_num,-1*cell2mat(all_SEM_81Q_ret_num),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Number of Cargoes');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);
% 
% figure('Name','Ant_ret_fractions','NumberTitle','off'), hold on,
% plot(box_sizes,all_mean_30Q_ant_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,all_mean_45Q_ant_frac,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,all_mean_65Q_ant_frac,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,all_mean_81Q_ant_frac,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,all_mean_30Q_ant_frac,cell2mat(all_SEM_30Q_ant_frac),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_45Q_ant_frac,cell2mat(all_SEM_45Q_ant_frac),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_65Q_ant_frac,cell2mat(all_SEM_65Q_ant_frac),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,all_mean_81Q_ant_frac,cell2mat(all_SEM_81Q_ant_frac),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% plot(box_sizes,-1*all_mean_30Q_ret_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor',colour_30Q);
% plot(box_sizes,-1*all_mean_45Q_ret_frac,'-','Color',colour_45Q,'Marker','o','MarkerEdgeColor',colour_45Q);
% plot(box_sizes,-1*all_mean_65Q_ret_frac,'-','Color',colour_65Q,'Marker','o','MarkerEdgeColor',colour_65Q);
% plot(box_sizes,-1*all_mean_81Q_ret_frac,'-','Color',colour_81Q,'Marker','o','MarkerEdgeColor',colour_81Q);
% errorbar(box_sizes,-1*all_mean_30Q_ret_frac,-1*cell2mat(all_SEM_30Q_ret_frac),'Color',colour_30Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_45Q_ret_frac,-1*cell2mat(all_SEM_45Q_ret_frac),'Color',colour_45Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_65Q_ret_frac,-1*cell2mat(all_SEM_65Q_ret_frac),'Color',colour_65Q,'LineWidth',2,'CapSize',0);
% errorbar(box_sizes,-1*all_mean_81Q_ret_frac,-1*cell2mat(all_SEM_81Q_ret_frac),'Color',colour_81Q,'LineWidth',2,'CapSize',0);
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Fraction of Cargoes');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% set(gca, 'FontSize',40);


%Previous Figure versions
% figure(), hold on,
% scatter(0.75-0.125+0.25*rand(numel(tot_ent_1um{1}),1),tot_ent_1um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(1.75-0.125+0.25*rand(numel(tot_ent_2um{1}),1),tot_ent_2um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(4.75-0.125+0.25*rand(numel(tot_ent_5um{1}),1),tot_ent_5um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(9.75-0.125+0.25*rand(numel(tot_ent_10um{1}),1),tot_ent_10um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(14.75-0.125+0.25*rand(numel(tot_ent_15um{1}),1),tot_ent_15um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([0.75,1.75,4.75,9.75,14.75],all_mean_30Q_tot_enters,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor','k');
% scatter(1.25-0.125+0.25*rand(numel(tot_ent_1um{4}),1),tot_ent_1um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(2.25-0.125+0.25*rand(numel(tot_ent_2um{4}),1),tot_ent_2um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(5.25-0.125+0.25*rand(numel(tot_ent_5um{4}),1),tot_ent_5um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(10.25-0.125+0.25*rand(numel(tot_ent_10um{4}),1),tot_ent_10um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(15.25-0.125+0.25*rand(numel(tot_ent_15um{4}),1),tot_ent_15um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([1.25,2.25,5.25,10.25,15.25],all_mean_81Q_tot_enters,'-','Color',colour_81Q,'Marker','square','MarkerEdgeColor','k');
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Number of Cargoes');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% 
% figure(), hold on,
% scatter(0.75-0.125+0.25*rand(numel(frac_ant_1um{1}),1),frac_ant_1um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(1.75-0.125+0.25*rand(numel(frac_ant_2um{1}),1),frac_ant_2um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(4.75-0.125+0.25*rand(numel(frac_ant_5um{1}),1),frac_ant_5um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(9.75-0.125+0.25*rand(numel(frac_ant_10um{1}),1),frac_ant_10um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(14.75-0.125+0.25*rand(numel(frac_ant_15um{1}),1),frac_ant_15um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([0.75,1.75,4.75,9.75,14.75],all_mean_30Q_ant_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor','k');
% scatter(1.25-0.125+0.25*rand(numel(frac_ant_1um{4}),1),frac_ant_1um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(2.25-0.125+0.25*rand(numel(frac_ant_2um{4}),1),frac_ant_2um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(5.25-0.125+0.25*rand(numel(frac_ant_5um{4}),1),frac_ant_5um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(10.25-0.125+0.25*rand(numel(frac_ant_10um{4}),1),frac_ant_10um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(15.25-0.125+0.25*rand(numel(frac_ant_15um{4}),1),frac_ant_15um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([1.25,2.25,5.25,10.25,15.25],all_mean_81Q_ant_frac,'-','Color',colour_81Q,'Marker','square','MarkerEdgeColor','k');
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Fraction Anterograde');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');
% 
% 
% figure(), hold on,
% scatter(0.75-0.125+0.25*rand(numel(frac_ret_1um{1}),1),frac_ret_1um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(1.75-0.125+0.25*rand(numel(frac_ret_2um{1}),1),frac_ret_2um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(4.75-0.125+0.25*rand(numel(frac_ret_5um{1}),1),frac_ret_5um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(9.75-0.125+0.25*rand(numel(frac_ret_10um{1}),1),frac_ret_10um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(14.75-0.125+0.25*rand(numel(frac_ret_15um{1}),1),frac_ret_15um{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([0.75,1.75,4.75,9.75,14.75],all_mean_30Q_ret_frac,'-','Color',colour_30Q,'Marker','o','MarkerEdgeColor','k');
% scatter(1.25-0.125+0.25*rand(numel(frac_ret_1um{4}),1),frac_ret_1um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(2.25-0.125+0.25*rand(numel(frac_ret_2um{4}),1),frac_ret_2um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(5.25-0.125+0.25*rand(numel(frac_ret_5um{4}),1),frac_ret_5um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(10.25-0.125+0.25*rand(numel(frac_ret_10um{4}),1),frac_ret_10um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% scatter(15.25-0.125+0.25*rand(numel(frac_ret_15um{4}),1),frac_ret_15um{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
% plot([1.25,2.25,5.25,10.25,15.25],all_mean_81Q_ret_frac,'-','Color',colour_81Q,'Marker','square','MarkerEdgeColor','k');
% xticks([1 2 5 10 15]);
% xlabel('Width (\mum)')
% ylabel('Fraction reterograde');
% publication_fig(0,0,1);
% box("off");
% pbaspect([2 1 1]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');