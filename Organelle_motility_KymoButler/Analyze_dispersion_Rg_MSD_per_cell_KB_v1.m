%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary (modified by Emily N. P. Prowse)
%Code #2 to run (after analyze trackmate motility).

%note for mac users, the pesky ._ files can cause this code to create empty
%cells, and you won't notice until the plotting code. To get rid of the ._
%files, go to the terminal and enter dot_clean directory/of/the/files/.

%Steps for running this code:
% 1.Update the addpaths (lines 11-16) and directories (lines 52-129) and variable DT.
% 2. Run the code for each k_choose by changing it manually. 


% close all; 
clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');


%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good

plus_proc_times=[];
average_fraction_of_outward_runs_per_cell=[];

%% Reset these parameters everytime you run the code !!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_MSD=[];
ab1=[];
kg=0;
variance_plus_traj=[];
variance_minus_traj=[];

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_MSD=[];
ab1=[];
kg=0;
variance_plus_traj=[];
variance_minus_traj=[];
proc_run_pos_2=[];
fl='*.mat';
an=0;
DT=0.12;
kf=0;
proc_runs=[];
proc_times=[];
rng=[];
tlp=-0.5;
avg_lr_plus=[];
avg_lr_minus=[];
fraction_proc_time=[]; %Emily added 20221121
fraction_plus_proc_time=[]; %Emily added
fraction_minus_proc_time=[]; %Emily added 20220927
fraction_diff_time=[]; %Emily added 20220927
fraction_stat_time=[]; %Emily added 20220927
fraction_plus_diff_time=[]; %Emily added 20221020
fraction_minus_diff_time=[]; %Emily added 20221020
fraction_proc=[]; %Emily added 20221121
fraction_proc_out=[]; %Emily added 20221121
fraction_proc_in=[]; %Emily added 20221121
fraction_diff=[]; %Emily added 20221121
fraction_diff_out=[]; %Emily added 20221121
fraction_diff_in=[]; %Emily added 20221121
fraction_stat=[]; %Emily added 20221121
fraction_proc_all=[]; %Emily added 20221121
fraction_proc_out_all=[]; %Emily added 20221121
fraction_proc_in_all=[]; %Emily added 20221121
fraction_diff_all=[]; %Emily added 20221121
fraction_diff_out_all=[]; %Emily added 20221121
fraction_diff_in_all=[]; %Emily added 20221121
fraction_stat_all=[]; %Emily added 20221121
plus_proc_times=[];
avg_lr_plus_per_cell=[];
avg_lr_minus_per_cell=[];
avg_frac_t_out_per_cell=[];
avg_frac_t_in_per_cell=[];
avg_frac_t_proc_per_cell=[];
avg_frac_t_diff_per_cell=[];
avg_frac_t_diff_out_per_cell=[]; %Emily added
avg_frac_t_diff_in_per_cell=[]; %Emily added
avg_frac_t_stat_per_cell=[]; %Emily added
proc_run=[]; %Emily added
diff_run=[]; %Emily added 20220927
stat_run=[]; %Emily added 20220927
proc_time=[]; %Emily added
stat_time=[]; %Emily added 20220927
diff_time=[]; % Emily added 20220927
gc=0; %Emily added 20230823

colour_index=1; %Emily added 20230206

tic
for k_choose = 4 %Don't forget to change min_lr based on your cargo! 1 for BDNF, 0.8 for mito
        
if k_choose == 1  % 30Q
    %BDNF
%     dat_dir = '/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_dir/';
    %BDNF+IFg
%     dat_dir = '/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_dir/';
   %Mito
%     dat_dir = '/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_dir/';
     %Lyso
%     dat_dir = '/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_dir/';
    %Lyso+IFg
    dat_dir = '/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_mats/';
    save_dir='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_pos/';
    save_dir_Rg='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_rg/';
    save_dir_MSD='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_msd/';
    save_dir_dirbias='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_dir/';

elseif k_choose == 2  % 45Q
%      %BDNF
%     dat_dir = '/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_dir/';
%Mito
%     dat_dir = '/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_dir/';
 %Lyso
    dat_dir = '/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_mats/';
    save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_pos/';
    save_dir_Rg='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_rg/';
    save_dir_MSD='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_msd/';
    save_dir_dirbias='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_dir/';
%     
elseif k_choose == 3  % 65Q
% BDNF
%     dat_dir = '/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_dir/';
%Mito
%     dat_dir = '/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_dir/';
% Lyso
    dat_dir = '/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_mats/';
    save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_pos/';
    save_dir_Rg='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_rg/';
    save_dir_MSD='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_msd/';
    save_dir_dirbias='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_dir/';
    
elseif k_choose == 4 % 81Q
    %BDNF
%     dat_dir = '/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_dir/';
%     %BDNF+IFg
%     dat_dir = '/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_dir/';
   %Mito
%     dat_dir = '/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_dir/';
    %Lyso
%     dat_dir = '/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/';
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_dir/';
%     %Lyso+IFg
    dat_dir = '/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_mats/';
    save_dir='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_pos/';
    save_dir_Rg='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_rg/';
    save_dir_MSD='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_msd/';
    save_dir_dirbias='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_dir/';
end
var=0;
fls=dir(fullfile(dat_dir,fl));
min_lr=1; %1 for BDNF, 0.8 for mito, let's approximate 1 for lyso also
numel(fls);
for k=1:numel(fls) %k is the current file
    pos=[];
    kg=kg+1;%still not sure why we need both kg and k, kg is the iteration number of the loop
    load(fullfile(dat_dir,fls(k).name)); %load filename
    ab = kbpos(~cellfun(@isempty,{kbpos.position})); % outputs a logical array of if the position is empty (0 for empty) or not.
    ab1{kg}=ab; %indexing by the iteration number of the loop to get the position array for that file
    rg_1=[];
    var=numel(ab)+var; %Tells me number of trajectories
    alphas=[]; %Emily added 20230720
    for j=1:numel(ab) %from 1 to the number of elements in the position file, which would be the number of trajectories for a file
        if numel(ab(j).position)>15 %only analyze trajectories that have more than 15 timepoints of positional data
        kf=kf+1; 
        curr_pos=ab(j).position;
        time=[1:length(curr_pos)]'.*DT; %Emily modified 20221019 so time and position have same length (changed numel to length)

            res=analyze_run_length_reversals_v5_per_cell_bias(time,curr_pos,min_lr,0,k_choose);
            proc_run=[proc_run;res.proc_run];
            proc_time=[proc_time;res.proc_time];
            avg_lr_plus=[avg_lr_plus;res.avg_lrplus];
            avg_lr_minus=[avg_lr_minus;res.avg_lrminus];
        
            fraction_proc_time=[fraction_proc_time;(sum(res.proc_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))];
            fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time(find(res.proc_run>0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
            fraction_minus_proc_time=[fraction_minus_proc_time;(sum(res.proc_time(find(res.proc_run<0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
            fraction_diff_time=[fraction_diff_time;(sum(res.diff_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find plus end runs
            fraction_stat_time=[fraction_stat_time;(sum(res.stat_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find plus end runs
            %Emily added this section 20221020 for testing whether min_lr is
            %set correctly
            fraction_plus_diff_time=[fraction_plus_diff_time;(sum(res.diff_time(find(res.diff_run>0)))/sum(res.diff_time))]; %Emily added 20221020
            fraction_minus_diff_time=[fraction_minus_diff_time;(sum(res.diff_time(find(res.diff_run<0)))/sum(res.diff_time))]; %Emily added 20221020
            
        rg_0=func_Rg_1D(curr_pos); %calculate radius of gyration
        rg_1=[rg_1,rg_0]; %adds the currently calculated rg to the previous one in a list. 
        dat_pos_tim{kf}=smooth(ab(j).position,10,'sgolay',1); %gives another smoothed position file for the current trajectory?
        elseif kf==0
        dat_pos_tim{1}=[];
        else %do not analyze trajectories with <15 timepoints of positions
        dat_pos_tim{kf}=[];
        end
    end %the end of the function analyzing trajectories
    pos=[pos,dat_pos_tim];
    position{k}=pos;

    save([save_dir, num2str(k_choose) 'Position_per_cell'],'position');
    
    % 1D MSD
    K=1:40; %time interval to calculate MSD
    Ndat=numel(pos); %Emily changed from numel(xk) to numel(pos) for this code
    res=msd_fun_v1_emmessedwithit(K,pos,DT,Ndat,0);% change the value to 0 to stop plotting
    pos=[]; %clearing the values stored in pos
    if res.slp>0 %Emily added 20230823 to get rid of cells with negative alpha
        gc=gc+1;
        res2{gc}=res; %saving this 1D MSD into res2, indexed by file number (cell)
    else
        disp(fls(k).name)
        disp('negative alpha')
    end
    rg_all1{k}=rg_1; %the list of rgs created on line 87 is saved into rg_all1, indexed also by file number (cell)
    clear res dat_pos_tim xk yk rg_2 Ndat;
    j=0;%returning the value of j to 0 so that the next trajectory will be analyzed from the first position
    
    avg_lr_plus_per_cell=[avg_lr_plus_per_cell;mean(avg_lr_plus(~isnan(avg_lr_plus)))]; %20230331
    avg_lr_minus_per_cell=[avg_lr_minus_per_cell;mean(avg_lr_minus(~isnan(avg_lr_minus)))]; %20230331

    fraction_proc=[fraction_proc;fraction_proc_time]; %Emily added 20221121
    fraction_proc_out=[fraction_proc_out;fraction_plus_proc_time]; %Emily added 20221121
    fraction_proc_in=[fraction_proc_in;fraction_minus_proc_time]; %Emily added 20221121
    fraction_diff=[fraction_diff;fraction_diff_time]; %Emily added 20221121
    fraction_diff_out=[fraction_diff_out;fraction_plus_diff_time]; %Emily added 20221121
    fraction_diff_in=[fraction_diff_in;fraction_minus_diff_time]; %Emily added 20221121
    fraction_stat=[fraction_stat;fraction_stat_time]; %Emily added 20221121
    
    avg_frac_t_proc_per_cell= [avg_frac_t_proc_per_cell;nanmean(fraction_proc_time)]; %Emily added
    avg_frac_t_out_per_cell= [avg_frac_t_out_per_cell;nanmean(fraction_plus_proc_time)]; %Emily added
    avg_frac_t_in_per_cell= [avg_frac_t_in_per_cell;nanmean(fraction_minus_proc_time)]; %Emily added
    avg_frac_t_diff_per_cell= [avg_frac_t_diff_per_cell;nanmean(fraction_diff_time)]; %Emily added
    avg_frac_t_diff_out_per_cell= [avg_frac_t_diff_out_per_cell;nanmean(fraction_plus_diff_time)]; %Emily added
    avg_frac_t_diff_in_per_cell= [avg_frac_t_diff_in_per_cell;nanmean(fraction_minus_diff_time)]; %Emily added
    avg_frac_t_stat_per_cell= [avg_frac_t_stat_per_cell;nanmean(fraction_stat_time)]; %Emily added
    
    
%     fraction_proc_per_cell{k}=fraction_proc; %Emily added 20221121
%     fraction_proc_out_per_cell{k}=fraction_proc_out; %Emily added 20221121
%     fraction_proc_in_per_cell{k}=fraction_proc_in; %Emily added 20221121
%     fraction_diff_per_cell{k}=fraction_diff; %Emily added 20221121
%     fraction_diff_out_per_cell{k}=fraction_diff_out; %Emily added 20221121
%     fraction_diff_in_per_cell{k}=fraction_diff_in; %Emily added 20221121
%     fraction_stat_per_cell{k}=fraction_stat; %Emily added 20221121
    avg_lr_plus=[]; %20230331
    avg_lr_minus=[]; %20230331
    fraction_proc_time=[]; %Emily added 20221121
    fraction_plus_proc_time=[]; %Emily added 20221121
    fraction_minus_proc_time=[]; %Emily added 20221121
    fraction_diff_time=[]; %Emily added 20221121
    fraction_plus_diff_time=[]; %Emily added 20221121
    fraction_minus_diff_time=[]; %Emily added 20221121
    fraction_stat_time=[]; %Emily added 20221121   
end

rg_1=[];

%Getting a list of the fractions processive, diffusive, and stationary for
%all trajectories
fraction_proc_all=[fraction_proc_all;fraction_proc]; %Emily added 20221121
fraction_proc_out_all=[fraction_proc_out_all;fraction_proc_out]; %Emily added 20221121
fraction_proc_in_all=[fraction_proc_in_all;fraction_proc_in]; %Emily added 20221121
fraction_diff_all=[fraction_diff_all;fraction_diff]; %Emily added 20221121
fraction_diff_out_all=[fraction_diff_out_all;fraction_diff_out]; %Emily added 20221121
fraction_diff_in_all=[fraction_diff_in_all;fraction_diff_in]; %Emily added 20221121
fraction_stat_all=[fraction_stat_all;fraction_stat]; %Emily added 20221121

fraction_proc=[]; %Emily added 20221121
fraction_proc_out=[]; %Emily added 20221121
fraction_proc_in=[]; %Emily added 20221121
fraction_diff=[]; %Emily added 20221121
fraction_diff_out=[]; %Emily added 20221121
fraction_diff_in=[]; %Emily added 20221121
fraction_stat=[]; %Emily added 20221121

% Defining variables to save
    rg_all=rg_all1;
    res1=res2;
    plus_rl_per_cell=avg_lr_plus_per_cell;
    minus_rl_per_cell=avg_lr_minus_per_cell;
    t_proc_per_cell=avg_frac_t_proc_per_cell;
    t_proc_out_per_cell=avg_frac_t_out_per_cell;
    t_proc_in_per_cell=avg_frac_t_in_per_cell;
    t_diff_per_cell=avg_frac_t_diff_per_cell;
    t_diff_out_per_cell=avg_frac_t_diff_out_per_cell;
    t_diff_in_per_cell= avg_frac_t_diff_in_per_cell;
    t_stat_per_cell=avg_frac_t_stat_per_cell;
    frac_proc_all=fraction_proc_all;
    frac_proc_out_all=fraction_proc_out_all;
    frac_proc_in_all=fraction_proc_in_all;
    frac_diff_all=fraction_diff_all;
    frac_diff_out_all=fraction_diff_out_all;
    frac_diff_in_all= fraction_diff_in_all;
    frac_stat_all=fraction_stat_all;
% Saving variables
    save([save_dir_Rg, 'Rg_per_cell'],'rg_all');
    save([save_dir_MSD, 'MSD_per_cell'],'res1');
    save([save_dir_dirbias, 'avg_rl_plus_per_cell'], 'plus_rl_per_cell');
    save([save_dir_dirbias, 'avg_rl_minus_per_cell'], 'minus_rl_per_cell');
    save([save_dir_dirbias, 'frac_proc_per_cell'], 't_proc_per_cell');
    save([save_dir_dirbias, 'frac_proc_out_per_cell'], 't_proc_out_per_cell');
    save([save_dir_dirbias, 'frac_proc_in_per_cell'], 't_proc_in_per_cell');
    save([save_dir_dirbias, 'frac_diff_per_cell'], 't_diff_per_cell');
    save([save_dir_dirbias, 'frac_diff_out_per_cell'], 't_diff_out_per_cell');
    save([save_dir_dirbias, 'frac_diff_in_per_cell'], 't_diff_in_per_cell');
    save([save_dir_dirbias, 'frac_stat_per_cell'], 't_stat_per_cell');
    save([save_dir_dirbias, 'proc_all'], 'frac_proc_all');
    save([save_dir_dirbias, 'proc_out_all'], 'frac_proc_out_all');
    save([save_dir_dirbias, 'proc_in_all'], 'frac_proc_in_all');
    save([save_dir_dirbias, 'diff_all'], 'frac_diff_all');
    save([save_dir_dirbias, 'diff_out_all'], 'frac_diff_out_all');
    save([save_dir_dirbias, 'diff_in_all'], 'frac_diff_in_all');
    save([save_dir_dirbias, 'stat_all'], 'frac_stat_all');
% Clearing all variables after running the code
    clear rg_all res1 res2 rg_all1;
    clear avg_lr_plus_per_cell avg_lr_minus_per_cell;
    clear fraction_proc_per_cell;
    clear fraction_proc_out_per_cell;
    clear fraction_proc_in_per_cell;
    clear fraction_diff_per_cell;
    clear fraction_diff_out_per_cell;
    clear fraction_diff_in_per_cell;
    clear fraction_stat_per_cell;
    
end

clear ab pos kf;

toc