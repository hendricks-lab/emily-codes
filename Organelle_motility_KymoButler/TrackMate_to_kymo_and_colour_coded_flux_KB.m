%% Emily Prowse 20240110
%This code turns trackmate trajectories into kymographs and calculates the
%flux of a in either direction within a specified distance box (eg. 5 um). 
% See the flux section to understand exactly how it's calculated.
% 
% This version colour codes the kymographs by flux and only works with
% sample images rather than all data.

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

% cmap = colormap(hot(100));
% colour_30Q=[0.5 0.5 0.5];
% colour_45Q=cmap(20,:);%[0.75 0 0];
% colour_65Q=cmap(25,:);%[0.5 0 0];
% colour_81Q=cmap(30,:);%[0.25 0 0];

cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_65Q=cmap(30,:);%[0.5 0 0];
cmap2= colormap(spring(100));
colour_45Q=[1.0000 0.4444 0.5556];%cmap2(50,:);%[0.75 0 0];
colour_81Q=cmap2(30,:);%[0.25 0 0];



count_stays_in_box=0;
count_ant_exits=0;
count_ret_exits=0;
count_entered_ant=0;
count_entered_ret=0;
count_did_not_enter=0;
count_started_in=0;
count_ent_then_stop=0;

% Mot_file='*.mat';
DT=0.12;    % single channel exposure time 120ms

%% Variables defined
for k_choose = 1

if k_choose == 1    % 30Q
    pth1=('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_mats/');

    pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_dir/';
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_dir/');

     Mot_file='1trackCoordinates_20230413_30Q_DIV8_lyso_s3_1002_kymo.mat'; %replaced with more representative kymograph 20240401
%      Mot_file='1trackCoordinates_20221104_30Q_DIV7_bdnf_qdot_s3_1003_kymo.mat';
%      Mot_file='1trackCoordinates_20230607_30Q_DIV7_mito_s1_1004_kymo.mat';
elseif k_choose == 2    % 45Q
%     pth1=('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_mats/'); 
%     pth1=('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_mats/');
    pth1=('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_mats/');


%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_Ifg_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_dir/');
    pth2=('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_dir/');

%     Mot_file='2trackCoordinates_20221118_45Q_DIV8_lyso_s2_1005_kymo.mat';
%     Mot_file='5trackCoordinates_20230717_30Q_DIV7_lyso_IFg_s1_1009_kymo.mat';
%     Mot_file='2trackCoordinates_20221216_45Q_BDNF_s1_1002_kymo.mat';
    Mot_file='5trackCoordinates_20230613_30Q_DIV7_IFg_bdnf_s2_1009_kymo.mat';
%     Mot_file='2trackCoordinates_20221118_45Q_DIV8_mito_s1_1007_kymo.mat';
elseif k_choose == 3   % 65Q
%     pth1=('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/');
    pth1=('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_mats/');
    
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/');
    pth2=('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_dir/');

%     Mot_file='3trackCoordinates_20210701_isoHD_65Q_DIV7_lyso_s2_1004_kymo.mat';
%     Mot_file='4trackCoordinates_20210705_isoHD_81Q_DIV7_lyso_s1_1012_kymo.mat';
    Mot_file='3trackCoordinates_20230621_65Q_DIV7_bdnf_s1_1003_kymo.mat';
%     Mot_file='4trackCoordinates_20221029_isoHD_81Q_DIV9_bdnf_qdot_s1_1010_kymo.mat';
%     Mot_file='3trackCoordinates_20230621_65Q_DIV7_mito_s3_1003_kymo.mat';
elseif k_choose == 4    % 81Q
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_mats/');
    pth1=('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_mats/');
%     pth1=('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_mats/');

%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_dir/');
    pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_dir/');
%     pth2=('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_dir/');

%      Mot_file='4trackCoordinates_20210705_isoHD_81Q_DIV7_lyso_s1_1012_kymo.mat';
%      Mot_file='6trackCoordinates_20230717_81Q_DIV7_lyso_IFg_s2_1004_kymo.mat';
%      Mot_file='4trackCoordinates_20221029_isoHD_81Q_DIV9_bdnf_qdot_s1_1010_kymo.mat';
     Mot_file='6trackCoordinates_20230624_81Q_DIV7_bdnf_IFg_s2_1007_kymo.mat';
%      Mot_file='4trackCoordinates_20230530_81Q_DIV6_mito_s1_1002_kymo.mat';
end

addpath(pth1);
addpath(pth2);

fl2='axon_length_per_cell.mat';

dat=[];

load(fullfile(pth2,fl2));
dmot=dir(fullfile(pth1,Mot_file));
ne_name=dmot.name(1:end-4);
lbls = filenames2labels(pth1);
idx=find(lbls == ne_name);
axon_length=axon_length{idx};

if numel(axon_length)~=numel(dmot)
    disp('You are missing files');
    keyboard
end

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
    Ndat=numel(kbpos);
    kp=0; 

for kd=1:Ndat
    kp=kp+1;
    colour=[0.5 0.5 0.5]; %resetting for colour coded trajectories

%     x_position=ab(kd).xk;%x position (um) %comment out for kymobutler
%     y_position=ab(kd).yk;%y position (um) %comment out for kymobutler
%     t_position=ab(kd).tk;% Time (s) %comment out for kymobutler
%     
% 
%     origin_x=ab(kd).origin_x; %comment out for kymobutler
%     origin_y=ab(kd).origin_y; %comment out for kymobutler
%     ax_l=axon_length{kd};

%create a 1D position coordinate from the x and y coordinates
%     X0 = [origin_x origin_y]; %cell center position [x y] %, comment out for kymobutler
%     r0{kd} = sqrt((x_position-X0(1)).^2 + (y_position-X0(2)).^2); % comment out for kymobutler

%flux analysis: how many trajectories cross a 2um box in the middle of the axon?
% position=r0{kd}; %without kymobutler
    t_position=kbpos(kd).rt;% Time (s)
    position=kbpos(kd).position; %with kymobutler
for it=1:5 
    box_width=2.5; %We want to look +/- box_width from the middle of the axon
    ax_seg=(axon_length*it)/5-(box_width+1); %segment of the axon
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

if count_ant_exits_tot>=1
    colour=[1 0 1];
elseif count_ret_exits_tot>=1
    colour=[0 1 0];
elseif count_stays_in_box_tot>=1
    colour=[0 0 1];
end

figure(k_choose*k*10)
    plot(position, t_position,'-','Color',colour,'linewidth',2)  
    xlabel('Position (\mum)'), ylabel('Time (s)'), hold on
    xlim([0,axon_length]);
    set(gca,'Ydir','reverse');
    set(gca,'FontSize',30);
    set(gca,'FontName','Arial');
    set(gca,'Box','off');

%     figure(k*1000)
%     plot(position{kd}, t_position{kd},'-','Color',cmap(kd,:),'linewidth',2)
%     xlabel('Position (\mum)'), ylabel('Time (s)'), hold on
%     xlim([0,axon_length]);
%     set(gca,'Ydir','reverse')

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

b1 = [(axon_length*1)/5-(box_width+1)-box_width (axon_length*1)/5-(box_width+1)+box_width];
b2 = [(axon_length*2)/5-(box_width+1)-box_width (axon_length*2)/5-(box_width+1)+box_width];
b3 = [(axon_length*3)/5-(box_width+1)-box_width (axon_length*3)/5-(box_width+1)+box_width];
b4 = [(axon_length*4)/5-(box_width+1)-box_width (axon_length*4)/5-(box_width+1)+box_width];
b5 = [(axon_length*5)/5-(box_width+1)-box_width (axon_length*5)/5-(box_width+1)+box_width];

figure(k_choose*k*10)
patch([b1(1) b1(1), b1(2) b1(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
patch([b2(1) b2(1), b2(2) b2(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
patch([b3(1) b3(1), b3(2) b3(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
patch([b4(1) b4(1), b4(2) b4(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
patch([b5(1) b5(1), b5(2) b5(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
pbaspect([2 1 1]);

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


clear x_position X0 y_position t_position_actual t_position

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

% save([save_dir, 'KB_total_enters_5um_newdef'],'tot_enters_per_cell_all_conds');
% save([save_dir, 'KB_frac_ant_ent_5um_newdef'],'frac_ant_enter_per_cell_all_conds');
% save([save_dir, 'KB_frac_ret_ent_5um_newdef'],'frac_ret_enter_per_cell_all_conds');
% save([save_dir, 'KB_start_in_5um_newdef'],'frac_starts_in_per_cell_all_conds');
% save([save_dir, 'KB_frac_ant_ex_5um_newdef'],'frac_ant_exits_per_cell_all_conds');
% save([save_dir, 'KB_frac_ret_ex_5um_newdef'],'frac_ret_exits_per_cell_all_conds');
% save([save_dir, 'KB_frac_stay_in_5um_newdef'],'frac_stays_in_per_cell_all_conds');

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

% fileprefix='20240110_isoHD_KB_bdnf_flux_5um_colour_coded_ex_';

% tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
% FolderName_1 = tempdir_1;   % Your destination folder
%Saving stats to a table
% stats_table=table(titles, rg_stats, alpha_stats,dirbias_stats, proc_stats, diff_stats, stat_stats);
% writetable(stats_table,fullfile(FolderName_1, [fileprefix,'.csv']),'WriteRowNames',true);

% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/';  % Your destination folder
% FolderName = tempdir;  % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%  FigHandle = FigList(iFig);
%  FigName  = get(FigHandle, 'Name');
%  savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end


