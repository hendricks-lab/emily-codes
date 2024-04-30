%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary

%Emily updated Sept 26/2022 to add the fraction of runs moving outward,
%inward, diffusively, or stationary per cell.
close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/raacampbell-notBoxPlot-7d90c27/code/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('E:\');
addpath('E:\New_General_Codes\');

%% Reset these parameters everytime you run the code !!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
k_choose = 4
    
if k_choose == 1    % 30Q 
    col1=[0.5 0.5 0.5];
      %BDNF
%     dat_dir = 'E:\Neuron\isoHD30Q\30Q_bdnf_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD30Q\30Q_bdnf_dir\';
     %Lyso
%    dat_dir = 'E:\Neuron\isoHD30Q\30Q_lyso_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD30Q\30Q_lyso_dir\';
    %Mito
   dat_dir = 'E:\Neuron\isoHD30Q\30Q_mito_mats\';
    save_dir_dirbias='E:\Neuron\isoHD30Q\30Q_mito_dir\';
    
    %BDNF
%     dat_dir = 'E:\Neurons\30Q_KB\30Q_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\30Q_KB\30Q_KB_bdnf_dir\';
    %Mito
%     dat_dir = 'E:\Neurons\30Q_KB\30Q_KB_mito_mats\';
%     save_dir_dirbias='E:\Neurons\30Q_KB\30Q_KB_mito_dir\';
     %Lyso
%     dat_dir = 'E:\Neurons\30Q_KB\30Q_KB_lyso_mats\'; %doesn't exist
%     save_dir_dirbias='E:\Neurons\30Q_KB\30Q_KB_lyso_dir\';
elseif k_choose == 2  % 45Q
    col1='b';
        %BDNF
    dat_dir = 'E:\Neuron\isoHD45Q\45Q_bdnf_mats\';
    save_dir_dirbias='E:\Neuron\isoHD45Q\45Q_bdnf_dir\';
     %Lyso
%    dat_dir = 'E:\Neuron\isoHD45Q\45Q_lyso_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD45Q\45Q_lyso_dir\';
%     Mito
%    dat_dir = 'E:\Neuron\isoHD45Q\45Q_mito_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD45Q\45Q_mito_dir\';

    %BDNF
%     dat_dir = 'E:\Neurons\45Q_KB\45Q_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\45Q_KB\45Q_KB_bdnf_dir\';
%     %BDNF+IFg
%     dat_dir = 'E:\Neurons\30Q_IFg_KB\30Q_IFg_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\30Q_IFg_KB\30Q_IFg_KB_bdnf_dir\';
    %Mito
%     dat_dir = 'E:\Neurons\45Q_KB\45Q_KB_mito_mats\';
%     save_dir_dirbias='E:\Neurons\45Q_KB\45Q_KB_mito_dir\';
   
elseif k_choose == 3    % 65Q
    col1=[0.5 0.5 0.5];
       %BDNF
    dat_dir = 'E:\Neuron\isoHD65Q\65Q_bdnf_mats\';
    save_dir_dirbias='E:\Neuron\isoHD65Q\65Q_bdnf_dir\';
%      %Lyso
%    dat_dir = 'E:\Neuron\isoHD65Q\65Q_lyso_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD65Q\65Q_lyso_dir\';
    %Mito
%    dat_dir = 'E:\Neuron\isoHD65Q\65Q_mito_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD65Q\65Q_mito_dir\';

        %BDNF
%     dat_dir = 'E:\Neurons\65Q_KB\65Q_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\65Q_KB\65Q_KB_bdnf_dir\';
    %81Q BDNF
%     dat_dir = 'E:\Neurons\81Q_KB\81Q_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\81Q_KB\81Q_KB_bdnf_dir\';
%     %Mito
%     dat_dir = 'E:\Neurons\65Q_KB\65Q_KB_mito_mats\';
%     save_dir_dirbias='E:\Neurons\65Q_KB\65Q_KB_mito_dir\';
    elseif k_choose == 4   % 81Q
    col1=[0.5 0.5 0.5];
    %     BDNF
%     dat_dir = 'E:\Neuron\isoHD81Q\81Q_bdnf_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD81Q\81Q_bdnf_dir\';
%      %Lyso
%    dat_dir = 'E:\Neuron\isoHD81Q\81Q_lyso_mats\';
%     save_dir_dirbias='E:\Neuron\isoHD81Q\81Q_lyso_dir\';
    %Mito
   dat_dir = 'E:\Neuron\isoHD81Q\81Q_mito_mats\';
    save_dir_dirbias='E:\Neuron\isoHD81Q\81Q_mito_dir\';

     %BDNF
%     dat_dir = 'E:\Neurons\81Q_KB\81Q_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\81Q_KB\81Q_KB_bdnf_dir\';
    %81Q IFg BDNF
%     dat_dir = 'E:\Neurons\81Q_IFg_KB\81Q_IFg_KB_bdnf_mats\';
%     save_dir_dirbias='E:\Neurons\81Q_IFg_KB\81Q_IFg_KB_bdnf_dir\';
    %Mito
%     dat_dir = 'E:\Neurons\81Q_KB\81Q_KB_mito_mats\';
%     save_dir_dirbias='E:\Neurons\81Q_KB\81Q_KB_mito_dir\';
end
ab1=[];
fl='*.mat';
DT=0.12;
all_proc_trajs=[];
all_diff_trajs=[];
proc_runs=[];
diff_runs=[];
stat_runs=[];
proc_alph1=[];
proc_alpha=[];
diff_alph1=[];
diff_alpha=[];
fraction_proc_time=[]; %Emily 20221102
fraction_plus_proc_time=[]; %Emily added
fraction_minus_proc_time=[]; %Emily added 20220927
fraction_diff_time=[]; %Emily added 20220927
fraction_stat_time=[]; %Emily added 20220927
fraction_plus_diff_time=[]; %Emily added 20221020
fraction_minus_diff_time=[]; %Emily added 20221020
proc_run=[]; %Emily added
delays = 1:20;

fls=dir(fullfile(dat_dir,fl));

min_lr_list = 0.2:0.2:3;
for iteration=1:numel(min_lr_list)%5%1:25 %Emily turned into a loop 20221031
    min_lr=min_lr_list(iteration);%0.1*iteration;
    kf=0; %Emily moved to here so that these reset for each iteration 20221116
    kg=0; %Emily moved to here so that these reset for each iteration 20221116
    
    proc_alpha_pc=[];
    kpr = 0;
    all_proc_runs_pc = [];
    diff_alpha_pc=[];
    kdr = 0;
    all_diff_runs_pc = [];
    for k=1:numel(fls) %k is the current file aka cell
        pos=[];
        kg=kg+1;%keeping track of the number of cells
        load(fullfile(dat_dir,fls(k).name), 'ab'); %load filename
        ab = ab(~cellfun(@isempty,{ab.position})); % outputs a logical array of if the position is empty (0 for empty) or not.
        ab1{k}=ab; %Emily replaced the indexing by kg with indexing by k instead 20221116
        rg_1=[];

        proc_alpha_pt = [];
        diff_alpha_pt = [];
        for j=1:numel(ab) %from 1 to the number of elements in the position file, which would be the number of trajectories for a file
            if numel(ab(j).position)>15 %only analyze trajectories that have more than 15 timepoints of positional data
                kf=kf+1; %keeping track of the number of trajectories for all cells

                curr_pos=ab(j).position;
                time=(1:length(curr_pos))'.*DT; %Emily modified 20221019 so time and position have same length (changed numel to length)
                res=analyze_run_length_reversals_v5_per_cell_bias(time,curr_pos,min_lr,0,k_choose);

                fraction_proc_time=[fraction_proc_time;sum(res.proc_time)/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))];
                fraction_diff_time=[fraction_diff_time;(sum(res.diff_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find fraction of diffusive runs
                fraction_stat_time=[fraction_stat_time;(sum(res.stat_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find fraction of stationary runs

                proc_runs=[proc_runs;res.proc_run]; %Emily added 20221116 to calculate number of runs
                diff_runs=[diff_runs;res.diff_run]; %Emily added 20221116 to calculate number of runs
                stat_runs=[stat_runs;res.stat_run]; %Emily added 20221116 to calculate number of runs

                fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time(find(res.proc_run>0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
                fraction_minus_proc_time=[fraction_minus_proc_time;(sum(res.proc_time(find(res.proc_run<0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs

                fraction_plus_diff_time=[fraction_plus_diff_time;(sum(res.diff_time(find(res.diff_run>0)))/sum(res.diff_time))]; %Emily added 20221020
                fraction_minus_diff_time=[fraction_minus_diff_time;(sum(res.diff_time(find(res.diff_run<0)))/sum(res.diff_time))]; %Emily added 20221020


                if numel(res.proc_traj)>=1 & ~isempty(res.proc_traj{1}),
                    [proc_msd_pt, proc_msd_std_pt, proc_msd_N_pt] = msd2(res.proc_traj,delays);
                    p_msd_pt=polyfit(log10(DT.*delays),log10(proc_msd_pt),1);
                    res.proc_alpha = p_msd_pt(1);

                    if ~isnan(res.proc_alpha)
                        proc_alpha_pt = [proc_alpha_pt; res.proc_alpha]; %Adam changed 20221115, for all the trajectories
                    end
                    
                    for jpr = 1:numel(res.proc_traj)
                    kpr=kpr+1;
                    all_proc_runs_pc{kpr}=res.proc_traj{jpr};
                    end
                end

                if numel(res.diff_traj)>=1 & ~isempty(res.diff_traj{1}),
                    [diff_msd_pt, diff_msd_std_pt, diff_msd_N_pt] = msd2(res.diff_traj,delays);
                    diff_p_msd_pt=polyfit(log10(DT.*delays),log10(diff_msd_pt),1);
                    res.diff_alpha = diff_p_msd_pt(1);

                    if ~isnan(res.diff_alpha)
                        diff_alpha_pt = [diff_alpha_pt; res.diff_alpha]; %Adam changed 20221115, for all the trajectories
                    end

                    for jdr = 1:numel(res.diff_traj)
                    kdr=kdr+1;
                    all_diff_runs_pc{kdr}=res.diff_traj{jdr};
                    end
                end

            else %do not analyze trajectories with <15 timepoints of positions.
            end
            j=0;%returning the value of j to 0 so that the next trajectory will be analyzed from the first position

        end %end of trajectory analysis

        proc_alpha_per_traj{k}=proc_alpha_pt;
        diff_alpha_per_traj{k}=diff_alpha_pt;
        fraction_proc_all{k}=fraction_proc_time;
        fraction_proc_out_all{k}=fraction_plus_proc_time;
        fraction_proc_in_all{k}=fraction_minus_proc_time;
        fraction_diff_all{k}=fraction_diff_time;
        fraction_diff_out_all{k}=fraction_plus_diff_time;
        fraction_diff_in_all{k}=fraction_minus_diff_time;
        fraction_stat_all{k}=fraction_stat_time;
        fraction_proc_all{k}=fraction_proc_time;
        fraction_proc_out_per_cell{k}=mean(fraction_plus_proc_time(find(~isnan(fraction_plus_proc_time))));
        fraction_proc_in_per_cell{k}=mean(fraction_minus_proc_time(find(~isnan(fraction_minus_proc_time))));
        fraction_diff_per_cell{k}=mean(fraction_diff_time(find(~isnan(fraction_diff_time))));
        fraction_diff_out_per_cell{k}=mean(fraction_plus_diff_time(find(~isnan(fraction_plus_diff_time))));
        fraction_diff_in_per_cell{k}=mean(fraction_minus_diff_time(find(~isnan(fraction_minus_diff_time))));
        fraction_stat_per_cell{k}=mean(fraction_stat_time(find(~isnan(fraction_stat_time))));
        num_proc_all{k}=numel(proc_runs); %Emily added 20221116
        num_diff_all{k}=numel(diff_runs); %Emily added 20221116
        num_stat_all{k}=numel(stat_runs); %Emily added 20221116


        fraction_proc_time=[];
        fraction_plus_proc_time=[];
        fraction_minus_proc_time=[];
        fraction_diff_time=[];
        fraction_plus_diff_time=[];
        fraction_minus_diff_time=[];
        fraction_stat_time=[];
        proc_run=[]; %Emily added 20221114
        proc_time=[]; %Emily added 20221114
        res=[]; %Emily added 20221114
        proc_runs=[]; %Emily added 20221116
        diff_runs=[]; %Emily added 20221116
        stat_runs=[]; %Emily added 20221116

        %calc. msd for all runs in cell
    [proc_msd_pc, proc_msd_std_pc, proc_msd_N_pc] = msd2(all_proc_runs_pc,delays);
    p_msd_pc=polyfit(log10(DT.*delays),log10(proc_msd_pc),1);
    results{iteration}.proc_alpha_per_cell{k} = p_msd_pc(1);

    [diff_msd_pc, diff_msd_std_pc, diff_msd_N_pc] = msd2(all_diff_runs_pc,delays);
    d_msd_pc=polyfit(log10(DT.*delays),log10(diff_msd_pc),1);
    results{iteration}.diff_alpha_per_cell{k} = d_msd_pc(1);

    end %End of cell analysis
    
    results{iteration}.proc_out=fraction_proc_out_all;
    results{iteration}.diff_out=fraction_diff_out_all;
    results{iteration}.proc_out_per_cell=fraction_proc_out_all;
    results{iteration}.diff_out_per_cell=fraction_diff_out_all;
    results{iteration}.proc_alpha_per_traj=proc_alpha_per_traj;
    results{iteration}.diff_alpha_per_traj=diff_alpha_per_traj; %Emily added 20221124
    results{iteration}.num_proc=num_proc_all;
    results{iteration}.num_diff=num_diff_all;
    results{iteration}.num_stat=num_stat_all;
    

    %     proc_alph_filtered=[];
    %     for cell=1:numel(proc_alpha_per_cell)
    %         if ~isempty(proc_alpha_per_cell(cell))
    %             proc_alph_filtered=[proc_alph_filtered;proc_alpha_per_cell(cell)];
    %         else
    %         end
    %     end
    %
    % %     diff_alph_filtered=[];
    % %     for cell=1:numel(diff_alpha_per_cell)
    % %         if ~isempty(diff_alpha_per_cell(cell))
    % %             diff_alph_filtered=[diff_alph_filtered;diff_alpha_per_cell(cell)];
    % %         else
    % %         end
    % %     end
    %
    %
        frac_diff_out_filtered=[];
        for cell=1:numel(fraction_diff_out_all)
            if ~isempty(fraction_diff_out_all)
                frac_diff_out_filtered=[frac_diff_out_filtered;fraction_diff_out_all(cell)];
            else
            end
        end
    
        frac_proc_out_filtered=[];
        for cell=1:numel(fraction_proc_out_all)
            if ~isempty(fraction_diff_out_all)
                frac_proc_out_filtered=[frac_proc_out_filtered;fraction_proc_out_all(cell)];
            else
            end
        end
        
        cdr=0;
        frac_diff_out_per_cell_filtered=[];
        for cell=1:numel(fraction_diff_out_per_cell)
            cdr=cdr+1;
            if ~isempty(fraction_diff_out_per_cell)
                frac_diff_out_per_cell_filtered=[frac_diff_out_per_cell_filtered;fraction_diff_out_per_cell(cdr)];
            else
            end
        end
        cpr=0;
        frac_proc_out_per_cell_filtered=[];
        for cell=1:numel(fraction_proc_out_per_cell)
            cpr=cpr+1;
            if ~isempty(fraction_diff_out_per_cell)
                frac_proc_out_per_cell_filtered=[frac_proc_out_per_cell_filtered;fraction_proc_out_per_cell(cpr)];
            else
            end
        end
    %
    %     proc_alpha_w_nans=cell2mat(proc_alph_filtered);
    %     %diff_alpha_w_nans=cell2mat(diff_alph_filtered);
    %     proc_alpha_mthd_2_w_nans=cell2mat(proc_alph_filtered);
    %     %diff_alpha_mthd_2_w_nans=cell2mat(diff_alph_filtered);
        frac_proc_out_mat_w_nans=cell2mat(frac_proc_out_filtered);
        frac_diff_out_mat_w_nans=cell2mat(frac_diff_out_filtered);
        frac_proc_out_per_cell_mat_w_nans=cell2mat(frac_proc_out_per_cell_filtered);
        frac_diff_out_per_cell_mat_w_nans=cell2mat(frac_diff_out_per_cell_filtered);
    %
    %     proc_alpha_mat=proc_alpha_w_nans(~isnan(proc_alpha_w_nans));
    %     %diff_alpha_mat=diff_alpha_w_nans(~isnan(diff_alpha_w_nans));
    %     proc_alpha_mthd_2_mat=proc_alpha_w_nans(~isnan(proc_alpha_mthd_2_w_nans));
    %     %diff_alpha_mthd_2_mat=diff_alpha_w_nans(~isnan(diff_alpha_mthd_2_w_nans));
        frac_proc_out_mat=frac_proc_out_mat_w_nans(~isnan(frac_proc_out_mat_w_nans));
        frac_diff_out_mat= frac_diff_out_mat_w_nans(~isnan(frac_diff_out_mat_w_nans));
        frac_proc_out_per_cell_mat=frac_proc_out_per_cell_mat_w_nans(~isnan(frac_proc_out_per_cell_mat_w_nans));
        frac_diff_out_per_cell_mat= frac_diff_out_per_cell_mat_w_nans(~isnan(frac_diff_out_per_cell_mat_w_nans));
    %
    %     avg_proc_alpha=mean(proc_alpha_mat);
    %     avg_diff_alpha=mean(diff_alpha_mat);
    %     avg_proc_alpha_mthd_2=mean(proc_alpha_mthd_2_mat);
    %     avg_diff_alpha_mthd_2=mean(diff_alpha_mthd_2_mat);
        avg_proc_out=mean(frac_proc_out_mat);
        avg_diff_out=mean(frac_diff_out_mat);
        avg_proc_out_per_cell=mean(frac_proc_out_per_cell_mat);
        avg_diff_out_per_cell=mean(frac_diff_out_per_cell_mat);
    %
    %     sem_proc_alpha=std(proc_alpha_mat)/sqrt(length(proc_alpha_mat));
    %     sem_diff_alpha=std(diff_alpha_mat)/sqrt(length(diff_alpha_mat));
    %     sem_proc_alpha_mthd_2=std(proc_alpha_mthd_2_mat)/sqrt(length(proc_alpha_mthd_2_mat));
    %     sem_diff_alpha_mthd_2=std(diff_alpha_mthd_2_mat)/sqrt(length(diff_alpha_mthd_2_mat));
        sem_frac_proc_out=std(frac_proc_out_mat)/sqrt(length(frac_proc_out_mat));
        sem_frac_diff_out=std(frac_diff_out_mat)/sqrt(length(frac_diff_out_mat));
    %
    %     %     nbins=50;
    %     %
    %     %     figure(200), hold on,
    %     %     histogram(proc_alpha1,nbins,'Normalization','probability')
    %     %     xlabel('Alpha for Processive Runs');
    %     %     ylabel('Frequency');
    %     %     publication_fig(0,0,1);
    %     %     nbins=50;
    %     %
    %     %     figure(200), hold on,
    %     %     histogram(diff_alpha1,nbins,'Normalization','probability')
    %     %     xlabel('Alpha for Processive Runs');
    %     %     ylabel('Frequency');
    %     %     publication_fig(0,0,1);
    %
    %     figure(50), hold on,
    %     plot(min_lr, avg_proc_alpha, 'ro');
    %     errorbar(min_lr,avg_proc_alpha,sem_proc_alpha,'k-');
    %     plot(min_lr,avg_diff_alpha,'bo');
    %     errorbar(min_lr,avg_diff_alpha,sem_diff_alpha,'k-');
    %     xlabel('Minimum Run Length (\mum)');
    %     ylabel('Average Run Alpha');
    %     publication_fig(0,0,1);
    %
    %
    %
    
        figure(10), hold on,
        plot(min_lr, avg_proc_out, 'ro');
        plot(min_lr,avg_diff_out,'b*');
        errorbar(min_lr,avg_proc_out,sem_frac_proc_out,'k-');
        errorbar(min_lr,avg_diff_out,sem_frac_diff_out,'k-');
        xlabel('Minimum Run Length (\mum)');
        ylabel('Fraction of Time of Runs Moving Outward');
        publication_fig(0,0,1);
        
        figure(11), hold on,
        scatter((min_lr-0.025)+0.05*rand(numel(frac_proc_out_mat),1),frac_proc_out_mat,25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3)
        scatter((min_lr-0.025)+0.05*rand(numel(frac_diff_out_mat),1),frac_diff_out_mat,25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3)
        plot(min_lr, avg_proc_out,'ko','markersize',10,'linewidth',2);
        plot(min_lr,avg_diff_out,'ks','markersize',10,'linewidth',2);
        xlabel('Minimum Run Length (\mum)');
        ylabel('Fraction of Time of Runs Moving Outward all runs');
        publication_fig(0,0,1);
        
        figure(12), hold on,
        scatter((min_lr-0.025)+0.05*rand(numel(frac_proc_out_per_cell_mat),1),frac_proc_out_per_cell_mat,25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3)
        scatter((min_lr-0.025)+0.05*rand(numel(frac_diff_out_per_cell_mat),1),frac_diff_out_per_cell_mat,25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3)
        plot(min_lr, avg_proc_out_per_cell,'ko','markersize',10,'linewidth',2);
        plot(min_lr,avg_diff_out_per_cell,'ks','markersize',10,'linewidth',2);
        xlabel('Minimum Run Length (\mum)');
        ylabel('Fraction of Time of Runs Moving Outward per cell');
        publication_fig(0,0,1);
        
%         figure(13), hold on,
%         plot(min_lr, avg_proc_out,'ko','markersize',10,'linewidth',2);
%         plot(min_lr,avg_diff_out,'ks','markersize',10,'linewidth',2);
%         scatter((min_lr-0.025)+0.05*rand(numel(frac_proc_out_mat),1),frac_proc_out_mat,25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3)
%         scatter((min_lr-0.025)+0.05*rand(numel(frac_diff_out_mat),1),frac_diff_out_mat,25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3)
%         xlabel('Minimum Run Length (\mum)');
%         ylabel('Fraction of Time of Runs Moving Outward per traj');
%         publication_fig(0,0,1);

        figure (101), hold on,
    proc_alpha_pt_all = [];
    for kc = 1:numel(results{iteration}.proc_alpha_per_traj)
        if ~isempty(results{iteration}.proc_alpha_per_traj{kc})
            scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.proc_alpha_per_traj{kc}),1),results{iteration}.proc_alpha_per_traj{kc},25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3);
            proc_alpha_pt_all = [proc_alpha_pt_all; results{iteration}.proc_alpha_per_traj{kc}];
            %plot(min_lr,mean(results{iteration}))
            %s.XJitterWidth = 0.5;
        end
    end
    plot(min_lr, mean(proc_alpha_pt_all),'ko','markersize',10,'linewidth',2)
    ylabel('alpha of processive runs (per traj.)')
    publication_fig(0,0,1);

    figure (102), hold on,
    diff_alpha_pt_all = [];
    for kc = 1:numel(results{iteration}.diff_alpha_per_traj)
        if ~isempty(results{iteration}.diff_alpha_per_traj{kc})
            scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.diff_alpha_per_traj{kc}),1),results{iteration}.diff_alpha_per_traj{kc},25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3);
            diff_alpha_pt_all = [diff_alpha_pt_all; results{iteration}.diff_alpha_per_traj{kc}];
            %plot(min_lr,mean(results{iteration}))
            %s.XJitterWidth = 0.5;
        end
    end
    plot(min_lr, mean(diff_alpha_pt_all),'ks','markersize',10,'linewidth',2)
    ylabel('alpha of diffusive runs (per traj.)')
    publication_fig(0,0,1);
    
    
%     figure (105), hold on,
%     proc_alpha_pt_all = [];
%     for kc = 1:numel(results{iteration}.proc_alpha_per_traj)
%         if ~isempty(results{iteration}.proc_alpha_per_traj{kc})
%             scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.proc_alpha_per_traj{kc}),1),results{iteration}.proc_alpha_per_traj{kc},25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3);
%             proc_alpha_pt_all = [proc_alpha_pt_all; results{iteration}.proc_alpha_per_traj{kc}];
%             %plot(min_lr,mean(results{iteration}))
%             %s.XJitterWidth = 0.5;
%         end
%     end
%     diff_alpha_pt_all = [];
%     for kc = 1:numel(results{iteration}.diff_alpha_per_traj)
%         if ~isempty(results{iteration}.diff_alpha_per_traj{kc})
%             scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.diff_alpha_per_traj{kc}),1),results{iteration}.diff_alpha_per_traj{kc},25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3);
%             diff_alpha_pt_all = [diff_alpha_pt_all; results{iteration}.diff_alpha_per_traj{kc}];
%             %plot(min_lr,mean(results{iteration}))
%             %s.XJitterWidth = 0.5;
%         end
%     end
%     plot(min_lr, mean(proc_alpha_pt_all),'ko','markersize',10,'linewidth',2)
%     ylabel('Alpha of Runs (per traj.)')
%     plot(min_lr, mean(diff_alpha_pt_all),'ks','markersize',10,'linewidth',2)
%     publication_fig(0,0,1);

    figure(103), hold on
    scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.proc_alpha_per_cell),1), cell2mat(results{iteration}.proc_alpha_per_cell),25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3)
    plot(min_lr, nanmean(cell2mat(results{iteration}.proc_alpha_per_cell)),'ko','markersize',10,'linewidth',2)
    ylabel('Alpha of Processive Runs (per cell)')
    publication_fig(0,0,1);

    figure(104), hold on
    scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.diff_alpha_per_cell),1), cell2mat(results{iteration}.diff_alpha_per_cell),25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3)
    plot(min_lr, nanmean(cell2mat(results{iteration}.diff_alpha_per_cell)),'ks','markersize',10,'linewidth',2)
    ylabel('Alpha of Diffusive Runs (per cell)')
    publication_fig(0,0,1);
    
    figure(106), hold on
    scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.proc_alpha_per_cell),1), cell2mat(results{iteration}.proc_alpha_per_cell),25,'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.3)
    scatter((min_lr-0.025)+0.05*rand(numel(results{iteration}.diff_alpha_per_cell),1), cell2mat(results{iteration}.diff_alpha_per_cell),25,'filled','MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3)
    plot(min_lr, mean(cell2mat(results{iteration}.proc_alpha_per_cell),'omitnan'),'ko','markersize',10,'linewidth',2)
    plot(min_lr, mean(cell2mat(results{iteration}.diff_alpha_per_cell),'omitnan'),'ks','markersize',10,'linewidth',2)
    publication_fig(0,0,1);
    ylabel('Alpha of Runs (per cell)')
    publication_fig(0,0,1);

    clear proc_alpha_per_traj proc_alph_filtered;
    clear diff_alpha_per_traj diff_alph_filtered;
    clear proc_alpha_mthd_2_per_cell proc_alph_mthd_2_filtered;
    clear diff_alpha_mthd_2_per_cell diff_alph_mthd_2_filtered;
    clear fraction_proc_per_cell fraction_proc_out_per_cell fraction_proc_in_per_cell_min_lrs fraction_diff_per_cell fraction_diff_out_per_cell fraction_diff_in_per_cell fraction_stat_per_cell;
    clear num_proc_per_cell num_diff_per_cell num_stat_per_cell;
end %end of min lr analysis
save([save_dir_dirbias,append(num2str(k_choose),'_mito_min_lr')], 'results');

toc
