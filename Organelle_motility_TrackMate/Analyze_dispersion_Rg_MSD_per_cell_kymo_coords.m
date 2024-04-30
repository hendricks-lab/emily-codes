%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary (modified by Emily N. P. Prowse)
%Code #2 to run (after analyze trackmate motility).

%note for mac users, the pesky ._ files can cause this code to create empty
%cells, and you won't notice until the plotting code. To get rid of the ._
%files, go to the terminal and enter dot_clean directory/of/the/files/.

%Steps for running this code:
% 1.Update the addpaths (lines 11-16) and directories (lines 52-129) and variable DT.
% 2. Run the code for each k_choose by changing it manually. 

%Feb 6/23 update, colour coded data based on experiment date, plots alpha & rg.

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

colour_index=1; %Emily added 20230206

tic
for k_choose = 2 %Don't forget to change min_lr based on your cargo! 0.5 for BDNF, 0.6 for mito, 0.8 for Lyso, also change K to 15 and change in 2d function for BDNF only
    
if k_choose == 1    % 18Q
    %Lyso-updated code not run for this condition
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_dir/';

% %BDNF-updated code not run for this condition
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_dir/';
    
elseif k_choose == 2  % 30Q
    %Lyso along kymo
    dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';

%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_dir/';
   %Lyso+IFg
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_dir/';
    %BDNF
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_dir/';

%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/';
    %BDNF+IFg
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_dir/';
   %Mito
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_dir/';

elseif k_choose == 3  % 45Q
%     %Lyso
%      dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_dir/';
     %Lyso along kymo
     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/';
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
%      %BDNF
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_dir/';
%Mito
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_dir/';
%     
elseif k_choose == 4  % 65Q
        %Lyso along kymo
    dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/';
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
%     %Lyso
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_dir/';
%BDNF
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_dir/';

    %Mito
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_dir/';
    
elseif k_choose == 5 % 81Q
     %     Lyso along kymo
    dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/';
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
    save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
% % % %     Lyso
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/';
    %Lyso+IFg
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_dir/';
    %BDNF
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_dir/';
%     %BDNF+IFg
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_dir/';
        %Mito
%     dat_dir = '/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_mats/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_dir/';
    
    

end
var=0;
fls=dir(fullfile(dat_dir,fl));
min_lr=0.5; %0.5 for BDNF, 0.6 for mito, 0.8 for lyso
numel(fls);
for k=1:numel(fls) %k is the current file
    pos=[];
    kg=kg+1;%still not sure why we need both kg and k, kg is the iteration number of the loop
    load(fullfile(dat_dir,fls(k).name)); %load filename
    ab = ab(~cellfun(@isempty,{ab.position})); % outputs a logical array of if the position is empty (0 for empty) or not.
    ab1{kg}=ab; %indexing by the iteration number of the loop to get the position array for that file
    rg_1=[];
    var=numel(ab)+var; %Tells me number of trajectories
    alphas=[]; %Emily added 20230720
    for j=1:numel(ab) %from 1 to the number of elements in the position file, which would be the number of trajectories for a cell
        if numel(ab(j).position)>15 %only analyze trajectories that have more than 15 timepoints of positional data
        kf=kf+1; %again not sure why we need this and j, now I sortof do because j is reset to zero at the end of the loop for the file while kf is not.
        xk{kf}=smooth(ab(j).xk,span,'sgolay',pwr); %smoothing the position function so we don't get too many sharp switches in direction which are artificial
        yk{kf}=smooth(ab(j).yk,span,'sgolay',pwr); %smoothing the position function so we don't get too many sharp switches in direction which are artificial
        %line 80 is calculating the 2D MSD, which is used for calculating
        %the per cell alpha value I'm pretty sure. Check the function if
        %you want to know how it is calculating the MSD.
        
        pos_on{kf}=ab(j).position;
        pos_off{kf}=ab(j).position_off;

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
            

        [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d_updated_em(xk{kf},yk{kf}, DT, k_choose);
        alphas=[alpha_1,alphas]; %Emily added 20230720
        if alpha_1>0 %not analyzing the tracks that are so flat they look
%         like they have negative slope I guess
        rg_0=func_Rg_Linda_v2(pos_on{kf},pos_off{kf}); 
        % (above) calculating the radius of gyration (rg) from the smoothed
        % position data of the trajectory
        rg_1=[rg_1,rg_0]; %adds the currently calculated rg to the previous one in a list. 
        dat_pos_tim{kf}=smooth(ab(j).position,10,'sgolay',1); %gives another smoothed position file for the current trajectory?
%         else %if alpha is less than or equal to zero, this condition will run

        % this condition seems to only be fulfilled when there are many
        % NaNs 
        
%         dat_pos_tim{kf}=[];
        end

        else %do not analyze trajectories with <15 timepoints of positions.
            dat_pos_tim{kf+1}=[]; %Emily added 20230817 when optimizing for the position projection along kymograph coordinates code
%             xk=[]; %Emily added 20230817 when optimizing for the position projection along kymograph coordinates code
%             yk=[]; %Emily added 20230817 when optimizing for the position projection along kymograph coordinates code
        end
    end %the end of the function analyzing trajectories
        alpha_2d{k}=[mean(alphas)]; %Emily added 20230720
        pos=[pos,dat_pos_tim];
        position{k}=pos(~cellfun('isempty',pos));
        axon_length_per_cell{k}=ab(1).axon_length; 

        save([save_dir, num2str(k_choose) 'Position_per_cell'],'position');

    
        % Calculate 1D MSD
        K=1:15; %Emily-Maybe we should use the 2D MSD instead of this 1D one, since I'm getting the NaNs in this analysis
        Ndat=numel(position{k}); %Emily-I think the problem is that xk is getting bigger with each run through the loop and therefore Ndat is increasing even though it shouldn't
        %below is a calculation of the 1D MSD
        res=msd_fun_v1_emmessedwithit(K,position{k},DT,Ndat,0);% change the value to 0 to stop plotting
        pos=[]; %clearing the values stored in pos
        res2{k}=res; %saving this 1D MSD into res2, indexed by file number (cell)
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
    
        avg_frac_t_proc_per_cell= [avg_frac_t_proc_per_cell;mean(fraction_proc_time, "omitnan")]; %Emily added
        avg_frac_t_out_per_cell= [avg_frac_t_out_per_cell;mean(fraction_plus_proc_time, "omitnan")]; %Emily added
        avg_frac_t_in_per_cell= [avg_frac_t_in_per_cell;mean(fraction_minus_proc_time, "omitnan")]; %Emily added
        avg_frac_t_diff_per_cell= [avg_frac_t_diff_per_cell;mean(fraction_diff_time, "omitnan")]; %Emily added
        avg_frac_t_diff_out_per_cell= [avg_frac_t_diff_out_per_cell;mean(fraction_plus_diff_time, "omitnan")]; %Emily added
        avg_frac_t_diff_in_per_cell= [avg_frac_t_diff_in_per_cell;mean(fraction_minus_diff_time, "omitnan")]; %Emily added
        avg_frac_t_stat_per_cell= [avg_frac_t_stat_per_cell;mean(fraction_stat_time, "omitnan")]; %Emily added
    
    
%       fraction_proc_per_cell{k}=fraction_proc; %Emily added 20221121
%       fraction_proc_out_per_cell{k}=fraction_proc_out; %Emily added 20221121
%       fraction_proc_in_per_cell{k}=fraction_proc_in; %Emily added 20221121
%       fraction_diff_per_cell{k}=fraction_diff; %Emily added 20221121
%       fraction_diff_out_per_cell{k}=fraction_diff_out; %Emily added 20221121
%       fraction_diff_in_per_cell{k}=fraction_diff_in; %Emily added 20221121
%       fraction_stat_per_cell{k}=fraction_stat; %Emily added 20221121
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

%Emily adding this section to colour data by experiment (20230206)
%     colour= [0 0.7 0.2; 0 0.5 0.5; 0.6 0 0.8; 0.8 0.5 0; 0.3 0 0.9];
% 
%     for n=1:numel(fls)
%         rg_cell=rg_all1{n};
%         
%         if isempty(rg_cell)
%         continue
%         else
%         [mean_rg,pci]=mle(rg_cell,'distribution','exp');
% 
%             if n==1
%  
%             figure(1), hold on
%             h1=notBoxPlot_nodotorline(mean_rg,k_choose+rand(1)*0.1,'style','line');
%             set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :), 'MarkerSize', 5, 'LineWidth', 2);
%             ylabel('Radius of Gyration (\mum)');
%             xlim([0 6]);
%             ylim([0 4]);
%         
%             figure(2), hold on
%             h2=notBoxPlot_nodotorline(res2{n}.slp,k_choose+rand(1)*0.1,'style','line');
%             set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :), 'MarkerSize', 5,'LineWidth', 2);
%             ylabel('\alpha');
%             xlim([0 6]);
%             ylim([0 4]);
% 
%             elseif fls(n).name(1:15)==fls(n-1).name(1:15)
%             figure(1), hold on
%             h1=notBoxPlot_nodotorline(mean_rg,k_choose+rand(1)*0.1,'style','line');
%             set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :), 'MarkerSize', 5, 'LineWidth', 2);
%             figure(2), hold on
%             h2=notBoxPlot_nodotorline(res2{n}.slp,k_choose+rand(1)*0.1,'style','line');
%             set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :), 'MarkerSize', 5,'LineWidth', 2);
%         else
%             colour_index=colour_index+1;
%             figure(1), hold on
%             h1=notBoxPlot_nodotorline(mean_rg,k_choose+rand(1)*0.1,'style','line');
%             set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :),'MarkerSize', 5,'LineWidth', 2);
%             figure(2), hold on
%             h2=notBoxPlot_nodotorline(res2{n}.slp,k_choose+rand(1)*0.1,'style','line');
%             set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour(colour_index, :), 'MarkerSize', 5,'LineWidth', 2);
%             end
%         end
%     end

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
    save([save_dir_MSD,'alphas'],'alpha_2d'); %Emily added 20230720
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
    save([save_dir_dirbias, 'axon_length_per_cell'],'axon_length_per_cell');
% Clearing all variables after running the code
    clear alpha_2d;%Emily added 20230720
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

clear ab xk yk pos kf pos_on pos_off;

toc

% figure()
% scatter(plus_rl_per_cell,cell2mat(axon_length_per_cell),'MarkerFaceColor', [1/k_choose 0 0]);
% figure()
% scatter(minus_rl_per_cell,cell2mat(axon_length_per_cell),'MarkerFaceColor', [1/k_choose 0 0]);