%% Emily Prowse 20230519
%This code turns trackmate trajectories into kymographs and calculates the
%flux of a in either direction within a 2 um box. See the flux section to
%understand exactly how it's calculated.

%Note: if you're planning to use the kymograph function, do not perform
%this for all trajectories-your computer will be mad. change the maximum
%kd to a maximum of 10 and k to a maximum of 5, only run one condition 
%(k_choose) at a time. You'll have to comment out the figures at the end
%and the variables that generate them.
clear all;
close all;
clc;

addpath('/Volumes/Emily_htt_2/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');

colour_30Q=[0.5 0.5 0.5];
colour_45Q=[0.75 0 0];
colour_65Q=[0.5 0 0];
colour_81Q=[0.25 0 0];

% colour_30Q=[0.5 0.5 0.5]; %30Q
% colour_45Q=[0.5 0.6 0.6]; %30Q IFg
% colour_65Q=[0.25 0 0]; %81Q
% colour_81Q=[0.25 0.1 0.1]; %81Q IFg

count_stays_in_box=0;
count_ant_exits=0;
count_ret_exits=0;
count_entered_ant=0;
count_entered_ret=0;
count_did_not_enter=0;
count_started_in=0;
count_ent_then_stop=0;

Mot_file='*.mat';
DT=0.12;    % single channel exposure time 120ms

% save_dir='/Volumes/Emily_htt_2/Neuron/mito_flux/';
save_dir='/Volumes/Emily_htt_2/Neuron/lyso_flux/';
% save_dir='/Volumes/Emily_htt_2/Neuron/bdnf_flux/';

%% Variables defined
for k_choose = 1:4

if k_choose == 1    % 30Q
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_mats/');
%     cd('/Volumes/Emily_2022/Fake_TM_flux_test_mat/');
elseif k_choose == 2    % 45Q
    cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_mats/');
%       cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_Ifg_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_mats/');
%     cd('/Volumes/Emily_2022/Fake_TM_flux_test_mat/');
elseif k_choose == 3   % 65Q
    cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_mats/');
%       cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats/');
%     cd('/Volumes/Emily_2022/Fake_TM_flux_test_mat/');
elseif k_choose == 4    % 81Q
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_mats/');
%     cd('/Volumes/Emily_2022/Fake_TM_flux_test_mat/');
end

dat=[];

dmot=dir(Mot_file); % Pick the trajectory files in .mat format

count_stays_in_box=0;
count_ant_exits=0;
count_ret_exits=0;
count_entered_ant=0;
count_entered_ret=0;
count_did_not_enter=0;
count_started_in=0;
count_ent_then_stop=0;
n=0;

for k=1:length(dmot) % k is a vector of length dmot
    display(k)
    load(dmot(k).name);%load everything in dmot
    [filepath,name,ext] = fileparts(dmot(k).name);
    display(dmot(k).name);  % display
    cmap = colormap(lines(100000));
    Ndat=numel(ab); %comment out for kymobutler
%     Ndat=numel(kbpos);
    kp=0; 

for kd=1:Ndat
    kp=kp+1;

    x_position=ab(kd).xk;%x position (um) %comment out for kymobutler
    y_position=ab(kd).yk;%y position (um) %comment out for kymobutler
    t_position=ab(kd).tk;% Time (s) %comment out for kymobutler
    

    origin_x=ab(kd).origin_x; %comment out for kymobutler
    origin_y=ab(kd).origin_y; %comment out for kymobutler
    axon_length=ab(kd).axon_length; %comment out for kymobutler

%create a 1D position coordinate from the x and y coordinates
    X0 = [origin_x origin_y]; %cell center position [x y] %, comment out for kymobutler
    r0{kd} = sqrt((x_position-X0(1)).^2 + (y_position-X0(2)).^2); % comment out for kymobutler
    %plotting kymograph from trackmate positons
%     time{kd}=t_position;k_ 
%     figure(k*10) %Plotting the kymograph
%     plot(r0{kd}, time{kd},'-','Color',cmap(kd,:),'linewidth',2)
%     xlabel('Position (\mum)'), ylabel('Time (s)'), hold on
%     xlim([0,axon_length]);
%     set(gca,'Ydir','reverse')

%flux analysis: how many trajectories cross a 2um box in the middle of the axon?
position=r0{kd}; %without kymobutler
% t_position=kbpos(kd).tk;% Time (s)
% position=kbpos(kd).position; %with kymobutler
for it=1:5
ax_seg=(axon_length*it)/5-2; %segment of the axon
box_width=0.5; %We want to look +/- box_width from the middle of the axon
    if any(position >= ax_seg-box_width & position <= ax_seg+box_width) %is it in the box?
        n=n+1;%count how many fulfill this condition
        idx_enters=min(find(position >= ax_seg-box_width & position <= ax_seg+box_width)); %what is the index of the earliest timepoint that the position enters the box
        if idx_enters==1
            count_started_in(it)=1;
            traj_starts_in_box=position(1)>= ax_seg-box_width & position(1) <=ax_seg+box_width; %does it start in the box?
            traj_enters_box_ret=0; %position(1) < mid_length-box_width %does it enter from the somal side?
            traj_enters_box_ant=0; %position(1) > mid_length+box_width %does it enter from the distal side?
        else
            traj_enters_box_ret=position(idx_enters-1) < ax_seg-box_width; %does it enter from the somal side?
            traj_enters_box_ant=position(idx_enters-1) > ax_seg+box_width; %does it enter from the distal side?
        end

        if traj_enters_box_ret ==1 %the trajectory actually entered from the retrograde side
            count_entered_ret(it)=1;
%             traj_exits_box=max(find(position(idx_enters:end)+ box_width> mid_length-box_width | position(idx_enters:end) < mid_length-box_width));
            traj_exits_box_ret=position(end) < ax_seg-box_width; %what is the index of the latest timepoint that the position exits the box
            traj_exits_box_ant=position(end) > ax_seg+box_width; %what is the index of the latest timepoint that the position exits the box
%             traj_stays_in_box=all(position(idx_enters:end)> mid_length-box_width & position(idx_enters:end)<mid_length+box_width);
            traj_ent_then_stop=all(position(end)>= ax_seg-box_width & position(end) <=ax_seg+box_width);
            if any(traj_exits_box_ret)
                count_ret_exits(it)=1;
            elseif any(traj_exits_box_ant)
                count_ant_exits(it)=1;
            elseif any(traj_ent_then_stop)
%                 count_ent_then_stop=count_ent_then_stop+1; %counting this as a separate case
                count_ant_exits(it)=1; %counting the ones that entered from the soma side as anterograde flux
            end

        elseif traj_enters_box_ant==1 %the trajectory entered the box from the anterograde side
            count_entered_ant(it)=1;
%             traj_exits_box=max(find(position(idx_enters:end)+ box_width> mid_length-box_width | position(idx_enters:end) < mid_length-box_width));
            traj_exits_box_ret=position(end) < ax_seg-box_width; %did it exit on the soma side?
            traj_exits_box_ant=position(end) > ax_seg+box_width; %did it exit on the distal side?
%             traj_stays_in_box=all(position(idx_enters:end)> mid_length-box_width & position(idx_enters:end)<mid_length+box_width);
            traj_ent_then_stop=all(position(end)>= ax_seg-box_width & position(end) <=ax_seg+box_width);
            if any(traj_exits_box_ret)
                count_ret_exits(it)=1;
            elseif any(traj_exits_box_ant)
                count_ant_exits(it)=1;
            elseif any(traj_ent_then_stop)
%                 count_ent_then_stop=count_ent_then_stop+1; %counting this as a separate case
                count_ret_exits(it)=1; %counting those that enter from the distal side as retrograde flux
            end
        elseif traj_starts_in_box==1
%             count_started_in=count_started_in+1; % was double counted
%             traj_exits_box=max(find(position(idx_enters:end)+ box_width> mid_length-box_width | position(idx_enters:end) < mid_length-box_width));
            traj_exits_box_ret=position(end) < ax_seg-box_width; %did it exit on the soma side?
            traj_exits_box_ant=position(end) > ax_seg+box_width; %did it exit on the distal side?
            traj_stays_in_box=all(position(end)>= ax_seg-box_width & position(end) <=ax_seg+box_width);
            if any(traj_exits_box_ret)
                count_ret_exits(it)=1;
            elseif any(traj_exits_box_ant)
                count_ant_exits(it)=1;
            elseif any(traj_stays_in_box)
                count_stays_in_box(it)=1;
            end

        end
    else
        count_did_not_enter(it)=1;
        clear traj_starts_in_box traj_enters_box_ret traj_enters_box_ant
    end
    %the below variables are not used, since it is rare to find a
    %trajectory in a given box, we add together all the segments in the
    %trajectory and reset the counts after each trajectory.
    
%     count_entered_ant_tot=sum(count_entered_ant)
%     total_enters_per_seg_calc(it)=(count_entered_ant+count_entered_ret+count_started_in);
%     frac_ant_exits_per_seg_calc(it)=count_ant_exits/total_enters_per_seg_calc(it);
%     frac_ret_exits_per_seg_calc(it)=count_ret_exits/total_enters_per_seg_calc(it);
%     frac_stays_in_per_seg_calc(it)=count_stays_in_box/total_enters_per_seg_calc(it);
%     frac_ant_enter_per_seg_calc(it)=count_entered_ant/total_enters_per_seg_calc(it);
%     frac_ret_enter_per_seg_calc(it)=count_entered_ret/total_enters_per_seg_calc(it);
%     frac_starts_in_per_seg_calc(it)=count_started_in/total_enters_per_seg_calc(it);
end
 
    count_entered_ret_tot=sum(count_entered_ret);
    count_entered_ant_tot=sum(count_entered_ant);
    count_started_in_tot=sum(count_started_in);
    count_ant_exits_tot=sum(count_ant_exits);
    count_ret_exits_tot=sum(count_ret_exits);
    count_stays_in_box_tot=sum(count_stays_in_box);

% tot_enters_per_traj(kd)=count_entered_ant_tot+count_entered_ret_tot+count_started_in_tot;
% tot_frac_ant_exits_per_traj(kd)=count_ant_exits_tot/tot_enters_per_traj(kd);
% tot_frac_ret_exits_per_traj(kd)=count_ret_exits_tot/tot_enters_per_traj(kd);
% tot_frac_stays_in_per_traj(kd)=count_stays_in_box_tot/tot_enters_per_traj(kd);
% tot_frac_ant_enter_per_traj(kd)=count_entered_ant_tot/tot_enters_per_traj(kd);
% tot_frac_ret_enter_per_traj(kd)=count_entered_ret_tot/tot_enters_per_traj(kd);
% tot_frac_starts_in_per_traj(kd)=count_started_in_tot/tot_enters_per_traj(kd);

tot_ant_exits_per_traj(kd)=count_ant_exits_tot;
tot_ret_exits_per_traj(kd)=count_ret_exits_tot;
tot_stays_in_per_traj(kd)=count_stays_in_box_tot;
tot_ant_enter_per_traj(kd)=count_entered_ant_tot;
tot_ret_enter_per_traj(kd)=count_entered_ret_tot;
tot_starts_in_per_traj(kd)=count_started_in_tot;

count_entered_ant=0; %reset these values so it actually counts each trajectory separately. Updated 20230814.
count_ret_exits=0;
count_ant_exits=0;
count_stays_in_box=0;
count_entered_ret=0;
count_started_in=0;
end

tot_enters_per_cell{k}=sum(tot_ant_enter_per_traj)+sum(tot_ret_enter_per_traj)+sum(tot_starts_in_per_traj);
tot_ant_exits_per_cell{k}=sum(tot_ant_exits_per_traj);
tot_ret_exits_per_cell{k}=sum(tot_ret_exits_per_traj);
tot_stays_in_per_cell{k}=sum(tot_stays_in_per_traj);
tot_ant_enter_per_cell{k}=sum(tot_ant_enter_per_traj);
tot_ret_enter_per_cell{k}=sum(tot_ret_enter_per_traj);
tot_starts_in_per_cell{k}=sum(tot_starts_in_per_traj);

frac_ant_exits_per_cell{k}=sum(tot_ant_exits_per_traj)./tot_enters_per_cell{k};
frac_ret_exits_per_cell{k}=sum(tot_ret_exits_per_traj)./tot_enters_per_cell{k};
frac_stays_in_per_cell{k}=sum(tot_stays_in_per_traj)./tot_enters_per_cell{k};
frac_ant_enter_per_cell{k}=sum(tot_ant_enter_per_traj)./tot_enters_per_cell{k};
frac_ret_enter_per_cell{k}=sum(tot_ret_enter_per_traj)./tot_enters_per_cell{k};
frac_starts_in_per_cell{k}=sum(tot_starts_in_per_traj)./tot_enters_per_cell{k};


clear axon_length x_position X0 y_position t_position_actual t_position

tot_ant_exits_per_traj=[];
tot_ret_exits_per_traj=[];
tot_stays_in_per_traj=[];
tot_ant_enter_per_traj=[];
tot_ret_enter_per_traj=[];
tot_starts_in_per_traj=[];

end

tot_enters_per_cell_mat=cell2mat(tot_enters_per_cell);
frac_ant_enter_per_cell_mat=cell2mat(frac_ant_enter_per_cell);
frac_ret_enter_per_cell_mat=cell2mat(frac_ret_enter_per_cell);
frac_starts_in_per_cell_mat=cell2mat(frac_starts_in_per_cell);
frac_ant_exits_per_cell_mat=cell2mat(frac_ant_exits_per_cell);
frac_ret_exits_per_cell_mat=cell2mat(frac_ret_exits_per_cell);
frac_stays_in_per_cell_mat=cell2mat(frac_stays_in_per_cell);

tot_enters_per_cell_value_mat=cell2mat(tot_enters_per_cell);

tot_enters_per_cell_all_conds{k_choose}=tot_enters_per_cell_mat;
frac_ant_enter_per_cell_all_conds{k_choose}=frac_ant_enter_per_cell_mat;
frac_ret_enter_per_cell_all_conds{k_choose}=frac_ret_enter_per_cell_mat;
frac_starts_in_per_cell_all_conds{k_choose}=frac_starts_in_per_cell_mat;
frac_ant_exits_per_cell_all_conds{k_choose}=frac_ant_exits_per_cell_mat;
frac_ret_exits_per_cell_all_conds{k_choose}=frac_ret_exits_per_cell_mat;
frac_stays_in_per_cell_all_conds{k_choose}=frac_stays_in_per_cell_mat;

total_enters_per_cell_value_all_conds{k_choose}=tot_enters_per_cell_value_mat;

mean_total_enters_per_cell_all_conds{k_choose}=mean(tot_enters_per_cell_mat,"omitnan");
mean_frac_ant_enter_per_cell_all_conds{k_choose}=mean(frac_ant_enter_per_cell_mat,"omitnan");
mean_frac_ret_enter_per_cell_all_conds{k_choose}=mean(frac_ret_enter_per_cell_mat,"omitnan");
mean_frac_starts_in_per_cell_all_conds{k_choose}=mean(frac_starts_in_per_cell_mat,"omitnan");
mean_frac_ant_exits_per_cell_all_conds{k_choose}=mean(frac_ant_exits_per_cell_mat,"omitnan");
mean_frac_ret_exits_per_cell_all_conds{k_choose}=mean(frac_ret_exits_per_cell_mat,"omitnan");
mean_frac_stays_in_per_cell_all_conds{k_choose}=mean(frac_stays_in_per_cell_mat,"omitnan");

avg_total_entered(k_choose)=mean(tot_enters_per_cell_mat,"omitnan");
avg_entered_ant(k_choose)=mean(frac_ant_enter_per_cell_mat,"omitnan");
avg_entered_ret(k_choose)=mean(frac_ret_enter_per_cell_mat,"omitnan");
avg_starts_in(k_choose)=mean(frac_starts_in_per_cell_mat,"omitnan");
avg_exits_ant(k_choose)=mean(frac_ant_exits_per_cell_mat,"omitnan");
avg_exits_ret(k_choose)=mean(frac_ret_exits_per_cell_mat,"omitnan");
avg_stays_in(k_choose)=mean(frac_stays_in_per_cell_mat,"omitnan");

sem_total_entered(k_choose)=std(tot_enters_per_cell_mat,"omitnan")/(sqrt(length(tot_enters_per_cell_mat)));
sem_entered_ant(k_choose)=std(frac_ant_enter_per_cell_mat,"omitnan")/(sqrt(length(frac_ant_enter_per_cell_mat)));
sem_entered_ret(k_choose)=std(frac_ret_enter_per_cell_mat,"omitnan")/(sqrt(length(frac_ret_enter_per_cell_mat)));
sem_starts_in(k_choose)=std(frac_starts_in_per_cell_mat,"omitnan")/(sqrt(length(frac_starts_in_per_cell_mat)));
sem_exits_ant(k_choose)=std(frac_ant_exits_per_cell_mat,"omitnan")/(sqrt(length(frac_ant_exits_per_cell_mat)));
sem_exits_ret(k_choose)=std(frac_ret_exits_per_cell_mat,"omitnan")/(sqrt(length(frac_ret_exits_per_cell_mat)));
sem_stays_in(k_choose)=std(frac_stays_in_per_cell_mat,"omitnan")/(sqrt(length(frac_stays_in_per_cell_mat)));

save([save_dir, 'total_enters_1um'],'tot_enters_per_cell_all_conds');
save([save_dir, 'frac_ant_ent_1um'],'frac_ant_enter_per_cell_all_conds');
save([save_dir, 'frac_ret_ent_1um'],'frac_ret_enter_per_cell_all_conds');
save([save_dir, 'start_in_1um'],'frac_starts_in_per_cell_all_conds');
save([save_dir, 'frac_ant_ex_1um'],'frac_ant_exits_per_cell_all_conds');
save([save_dir, 'frac_ret_ex_1um'],'frac_ret_exits_per_cell_all_conds');
save([save_dir, 'frac_stay_in_1um'],'frac_stays_in_per_cell_all_conds');

tot_enters_per_cell=[];
tot_ant_exits_per_cell=[];
tot_ret_exits_per_cell=[];
tot_stays_in_per_cell=[];
tot_ant_enter_per_cell=[];
tot_ret_enter_per_cell=[];
tot_starts_in_per_cell=[];

frac_ant_exits_per_cell=[];
frac_ret_exits_per_cell=[];
frac_stays_in_per_cell=[];
frac_ant_enter_per_cell=[];
frac_ret_enter_per_cell=[];
frac_starts_in_per_cell=[];

end

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
xticklabels({'Anterograde flux', 'Retrograde flux','Stationary'});
ylabel('Fraction of Cargoes that enter the box');
publication_fig(0,0,1);

figure('Name','ant_ret_stat_entry_bar','NumberTitle','off'), hold on, %4 conditions
enter_fractions=[avg_entered_ant(1) avg_entered_ant(2) avg_entered_ant(3) avg_entered_ant(4); avg_entered_ret(1) avg_entered_ret(2) avg_entered_ret(3) avg_entered_ret(4); avg_starts_in(1) avg_starts_in(2) avg_starts_in(3) avg_starts_in(4)];
enter_err_bars=[sem_entered_ant(1) sem_entered_ant(2) sem_entered_ant(3) sem_entered_ant(4); sem_entered_ret(1) sem_entered_ret(2) sem_entered_ret(3) sem_entered_ret(4); sem_starts_in(1) sem_starts_in(2) sem_starts_in(3) sem_starts_in(4)];
b2=bar(enter_fractions,'grouped');
b2(1).FaceColor=colour_30Q;
b2(2).FaceColor=colour_45Q;
b2(3).FaceColor=colour_65Q;
b2(4).FaceColor=colour_81Q;
b2(1).FaceAlpha=0.5;
b2(2).FaceAlpha=0.5;
b2(3).FaceAlpha=0.5;
b2(4).FaceAlpha=0.5;
hold on
[ngroups,nbars]= size(enter_fractions);
groupwidth= min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
   % Calculate center of each bar
   x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
   errorbar(x, enter_fractions (:,i), enter_err_bars (:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
xticks([1 2 3]);
xticklabels({'Anterograde flux', 'Retrograde flux','Stationary'});
ylabel('Fraction of Cargoes that enter the box');
publication_fig(0,0,1);
xticks([1 2 3]);
xticklabels({'Somal Entry', 'Distal Entry','Started In'});
ylabel('Fraction of Cargoes that enter the box');
publication_fig(0,0,1);

figure('Name','ant_vs_ret_entry','NumberTitle','off'), hold on, %4 conditions
ant_vs_ret_fractions=[avg_entered_ant(1)/avg_entered_ret(1) avg_entered_ant(2)/avg_entered_ret(2) avg_entered_ant(3)/avg_entered_ret(3) avg_entered_ant(4)/avg_entered_ret(4)];
b3=bar(ant_vs_ret_fractions);
xticklabels({'30Q','45Q', '65Q', '81Q'});
% xticklabels({'30Q','30Q IFγ', '81Q', '81Q IFγ'});
ylabel('Anterograde/Retrograde Flux Entry');
publication_fig(0,0,1);

figure('Name','ant_vs_ret_exit','NumberTitle','off'), hold on, %4 conditions
ant_vs_ret_fractions=[avg_exits_ant(1)/avg_exits_ret(1) avg_exits_ant(2)/avg_exits_ret(2) avg_exits_ant(3)/avg_exits_ret(3) avg_exits_ant(4)/avg_exits_ret(4)];
b3=bar(ant_vs_ret_fractions);
xticklabels({'30Q','45Q', '65Q', '81Q'});
% xticklabels({'30Q','30Q IFγ', '81Q', '81Q IFγ'});
ylabel('Anterograde/Retrograde Flux Exit');
publication_fig(0,0,1);

figure('Name','total_entered','NumberTitle','off'), hold on, %4 conditions
scatter(1-0.125+0.25*rand(numel(total_enters_per_cell_value_all_conds{1}),1),total_enters_per_cell_value_all_conds{1},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2), hold on;
scatter(2-0.125+0.25*rand(numel(total_enters_per_cell_value_all_conds{2}),1),total_enters_per_cell_value_all_conds{2},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(total_enters_per_cell_value_all_conds{3}),1),total_enters_per_cell_value_all_conds{3},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(total_enters_per_cell_value_all_conds{4}),1),total_enters_per_cell_value_all_conds{4},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
plot(1,mean(total_enters_per_cell_value_all_conds{1},"omitnan"),'k_','MarkerSize',30);
plot(2,mean(total_enters_per_cell_value_all_conds{2},"omitnan"),'k_','MarkerSize',30);
plot(3,mean(total_enters_per_cell_value_all_conds{3},"omitnan"),'k_','MarkerSize',30);
plot(4,mean(total_enters_per_cell_value_all_conds{4},"omitnan"),'k_','MarkerSize',30);
xlim([0 5]);
xticks([1 2 3 4]);
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
xticklabels({'30Q','45Q', '65Q', '81Q'});
% xticklabels({'30Q','30Q IFγ', '81Q', '81Q IFγ'});
ylabel('Total Entered per cell');
publication_fig(0,0,1);

exit_all{1}=frac_ant_exits_per_cell_all_conds{1};
exit_all{2}=frac_ant_exits_per_cell_all_conds{2};
exit_all{3}=frac_ant_exits_per_cell_all_conds{3};
exit_all{4}=frac_ant_exits_per_cell_all_conds{4};
exit_all{5}=frac_ret_exits_per_cell_all_conds{1};
exit_all{6}=frac_ret_exits_per_cell_all_conds{2};
exit_all{7}=frac_ret_exits_per_cell_all_conds{3};
exit_all{8}=frac_ret_exits_per_cell_all_conds{4};
exit_all{9}=frac_stays_in_per_cell_all_conds{1};
exit_all{10}=frac_stays_in_per_cell_all_conds{2};
exit_all{11}=frac_stays_in_per_cell_all_conds{3};
exit_all{12}=frac_stays_in_per_cell_all_conds{4};

figure('Name','ant_ret_stat_exit_scatter','NumberTitle','off'), hold on,
% violin(exit_all,'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q;colour_30Q;colour_45Q;colour_65Q;colour_81Q;colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
scatter(1-0.125+0.25*rand(numel(exit_all{1}),1),exit_all{1},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(exit_all{2}),1),exit_all{2},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(exit_all{3}),1),exit_all{3},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(exit_all{4}),1),exit_all{4},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(exit_all{5}),1),exit_all{5},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(exit_all{6}),1),exit_all{6},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(exit_all{7}),1),exit_all{7},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(exit_all{8}),1),exit_all{8},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
scatter(11-0.125+0.25*rand(numel(exit_all{9}),1),exit_all{9},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(12-0.125+0.25*rand(numel(exit_all{10}),1),exit_all{10},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(13-0.125+0.25*rand(numel(exit_all{11}),1),exit_all{11},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(14-0.125+0.25*rand(numel(exit_all{12}),1),exit_all{12},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
plot(1,mean(exit_all{1},"omitnan"),'k_','MarkerSize',20);
plot(2,mean(exit_all{2},"omitnan"),'k_','MarkerSize',20);
plot(3,mean(exit_all{3},"omitnan"),'k_','MarkerSize',20);
plot(4,mean(exit_all{4},"omitnan"),'k_','MarkerSize',20);
plot(6,mean(exit_all{5},"omitnan"),'k_','MarkerSize',20);
plot(7,mean(exit_all{6},"omitnan"),'k_','MarkerSize',20);
plot(8,mean(exit_all{7},"omitnan"),'k_','MarkerSize',20);
plot(9,mean(exit_all{8},"omitnan"),'k_','MarkerSize',20);
plot(11,mean(exit_all{9},"omitnan"),'k_','MarkerSize',20);
plot(12,mean(exit_all{10},"omitnan"),'k_','MarkerSize',20);
plot(13,mean(exit_all{11},"omitnan"),'k_','MarkerSize',20);
plot(14,mean(exit_all{12},"omitnan"),'k_','MarkerSize',20);
xlim([0 15]);
ylabel('Fraction that enter the box');
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');


enter_all{1}=frac_ant_enter_per_cell_all_conds{1};
enter_all{2}=frac_ant_enter_per_cell_all_conds{2};
enter_all{3}=frac_ant_enter_per_cell_all_conds{3};
enter_all{4}=frac_ant_enter_per_cell_all_conds{4};
enter_all{5}=frac_ret_enter_per_cell_all_conds{1};
enter_all{6}=frac_ret_enter_per_cell_all_conds{2};
enter_all{7}=frac_ret_enter_per_cell_all_conds{3};
enter_all{8}=frac_ret_enter_per_cell_all_conds{4};
enter_all{9}=frac_stays_in_per_cell_all_conds{1};
enter_all{10}=frac_stays_in_per_cell_all_conds{2};
enter_all{11}=frac_stays_in_per_cell_all_conds{3};
enter_all{12}=frac_stays_in_per_cell_all_conds{4};

figure('Name','ant_ret_stat_entry_scatter','NumberTitle','off'), hold on,
% violin(enter_all,'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q;colour_30Q;colour_45Q;colour_65Q;colour_81Q;colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
scatter(1-0.125+0.25*rand(numel(enter_all{1}),1),enter_all{1},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(enter_all{2}),1),enter_all{2},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(enter_all{3}),1),enter_all{3},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(enter_all{4}),1),enter_all{4},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
scatter(6-0.125+0.25*rand(numel(enter_all{5}),1),enter_all{5},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(7-0.125+0.25*rand(numel(enter_all{6}),1),enter_all{6},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(8-0.125+0.25*rand(numel(enter_all{7}),1),enter_all{7},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(9-0.125+0.25*rand(numel(enter_all{8}),1),enter_all{8},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
scatter(11-0.125+0.25*rand(numel(enter_all{9}),1),enter_all{9},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2);
scatter(12-0.125+0.25*rand(numel(enter_all{10}),1),enter_all{10},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(13-0.125+0.25*rand(numel(enter_all{11}),1),enter_all{11},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(14-0.125+0.25*rand(numel(enter_all{12}),1),enter_all{12},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
plot(1,mean(enter_all{1},"omitnan"),'k_','MarkerSize',20);
plot(2,mean(enter_all{2},"omitnan"),'k_','MarkerSize',20);
plot(3,mean(enter_all{3},"omitnan"),'k_','MarkerSize',20);
plot(4,mean(enter_all{4},"omitnan"),'k_','MarkerSize',20);
plot(6,mean(enter_all{5},"omitnan"),'k_','MarkerSize',20);
plot(7,mean(enter_all{6},"omitnan"),'k_','MarkerSize',20);
plot(8,mean(enter_all{7},"omitnan"),'k_','MarkerSize',20);
plot(9,mean(enter_all{8},"omitnan"),'k_','MarkerSize',20);
plot(11,mean(enter_all{9},"omitnan"),'k_','MarkerSize',20);
plot(12,mean(enter_all{10},"omitnan"),'k_','MarkerSize',20);
plot(13,mean(enter_all{11},"omitnan"),'k_','MarkerSize',20);
plot(14,mean(enter_all{12},"omitnan"),'k_','MarkerSize',20);
xlim([0 15]);
ylabel('Fraction that enter the box');
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

ant_vs_ret_all{1}=frac_ant_exits_per_cell_all_conds{1}./frac_ret_exits_per_cell_all_conds{1};
ant_vs_ret_all{2}=frac_ant_exits_per_cell_all_conds{2}./frac_ret_exits_per_cell_all_conds{2};
ant_vs_ret_all{3}=frac_ant_exits_per_cell_all_conds{3}./frac_ret_exits_per_cell_all_conds{3};
ant_vs_ret_all{4}=frac_ant_exits_per_cell_all_conds{4}./frac_ret_exits_per_cell_all_conds{4};

figure('Name','ant_vs_ret_scatter','NumberTitle','off'), hold on,
scatter(1-0.125+0.25*rand(numel(ant_vs_ret_all{1}),1),ant_vs_ret_all{1},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_30Q,'LineWidth',2), hold on;
scatter(2-0.125+0.25*rand(numel(ant_vs_ret_all{2}),1),ant_vs_ret_all{2},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_45Q,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(ant_vs_ret_all{3}),1),ant_vs_ret_all{3},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_65Q,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(ant_vs_ret_all{4}),1),ant_vs_ret_all{4},25,'MarkerFaceColor','none','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',colour_81Q,'LineWidth',2);
plot(1,mean(ant_vs_ret_all{1}(~isinf(ant_vs_ret_all{1})), "omitnan"),'k_','MarkerSize',30);
plot(2,mean(ant_vs_ret_all{2}(~isinf(ant_vs_ret_all{2})), "omitnan"),'k_','MarkerSize',30);
plot(3,mean(ant_vs_ret_all{3}(~isinf(ant_vs_ret_all{3})), "omitnan"),'k_','MarkerSize',30);
plot(4,mean(ant_vs_ret_all{4}(~isinf(ant_vs_ret_all{4})), "omitnan"),'k_','MarkerSize',30);
ylabel('Fraction anterograde/retrograde');
xlim([0 5]);
xticks([1 2 3 4]);
xticklabels({'30Q','45Q', '65Q', '81Q'});
% xticklabels({'30Q','30Q IFγ', '81Q', '81Q IFγ'});
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

fileprefix='20230904_isoHD_lyso_flux_1um_box';
% 
% tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
% FolderName_1 = tempdir_1;   % Your destination folder
%Saving stats to a table
% stats_table=table(titles, rg_stats, alpha_stats,dirbias_stats, proc_stats, diff_stats, stat_stats);
% writetable(stats_table,fullfile(FolderName_1, [fileprefix,'.csv']),'WriteRowNames',true);

tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/';  % Your destination folder
FolderName = tempdir;  % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
 FigHandle = FigList(iFig);
 FigName  = get(FigHandle, 'Name');
 savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
 saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end


