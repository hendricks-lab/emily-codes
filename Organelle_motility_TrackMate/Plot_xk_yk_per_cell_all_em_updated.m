%% Plot x,y per cell with the center:
%% Abdullah R. Chaudhary - 05/16/2020

%Emily updated March 2022

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Volumes/Emily_htt/New_General_codes/');
addpath('/Volumes/Emily_htt/AC_codes_epmodified_20200603/');
addpath('/Volumes/Emily_htt/New_General_codes/plotSpread/');

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good

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

tic
for k_choose = 2
    
if k_choose == 1    % WT - Rab5
    col1=[0.5 0.5 0.5];
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_mats/');
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_mats/');
    fl='*.mat';
%     fl='1Spots in tracks statistics_20210123_WT_lyso_s1_1004.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    nm='WT';
    
elseif k_choose == 2   % S421D - Rab5
    col1=[0 0.4 0.4];
    cd('/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_mats/');
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_mats/');
    fl='*.mat';
    fl='2Spots in tracks statistics_20201025_S421D_r5_s4_1008.mat';
%     fl='2Spots in tracks statistics_20210116_S421D_lyso_s3_1009.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=0.5;
    nm='S421D';
    
end

%save_fig_final='/Volumes/Emily_htt/';

fls=dir(fl);
pos=[];  
for k=1:numel(fls)
% for k=1:30 %k is the cells/files
    kg=kg+1;
    %this section is because of an error in which ._ is added to the beginning of the first filename
    fln=fls(k).name;
    if fln(1:2)=="._"
    fln=fln(3:end);
    display(fln)
    load(fln);
    else
    load(fls(k).name);
    display(fls(k).name);
    end
    ab = ab(~cellfun(@isempty,{ab.position}));
    ab1{kg}=ab;
    rg_1=[];
    
%     for j=1:numel(ab)
    for j=1:100 %j is the trajectory number
%         numel(ab(j).position)
        if numel(ab(j).position)>15 %Emily changed this value from 75 to 0 on 15/08/2020
        kf=kf+1;
        xk{kf}=ab(j).xk;
        yk{kf}=ab(j).yk;
%         xk{kf}=ab(j).xk-(ab(j).xk(1)+randi(30,1,1)); %Emily modifying for plotting trajectories
%         yk{kf}=ab(j).yk-(ab(j).yk(1)+randi(10,1,1));
%         pos_tim{kf}=ab(j).position+randi(40,1,1); %Idea to have constant spacing between trajectories (base it on kf)
%         tim{kf}=ab(j).tk+randi(100,1,1);

        %Emily removing to try something else
%         if kf<26
%         pos_tim{kf}=ab(j).position+2*kf; %Idea to have constant spacing between trajectories (base it on kf)
%         tim{kf}=ab(j).tk+randi(120,1,1);%Added randi150 for lyso WT, 200 for lyso S421D, 120 for both early conditions
%         else
%         end
        %Emily uncommented
        pos_tim{kf}=ab(j).position;
        tim{kf}=ab(j).tk;
        xk{kf}=smooth(ab(j).xk,span,'sgolay',pwr); %Emily removed
%         smoothing on 24/03/2022
        yk{kf}=smooth(ab(j).yk,span,'sgolay',pwr);
        
        [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d (xk{kf},yk{kf}, DT, k_choose);
        
        if alpha_1>0 
        rg_0=func_Rg_Linda_v2(xk{kf},yk{kf});
        rg_1=[rg_1,rg_0];
        dat_pos_tim{kf}=smooth(ab(j).position,30,'sgolay',1);
        xk_yk{kf}=[xk{kf},yk{kf}];
        else
        rg_1=[];
        dat_pos_tim{kf}=[];
        xk_yk{kf}=[];
        end
        
        else
        end
      
    end
    pos=[pos,dat_pos_tim];
    position{k}=pos;
    
    %figure(k), hold on, 
    %cellfun(@(x,y) plot(x,y, 'Color','k'), xk_yk);
%     fig=figure(k); %Em commented out to try another style of figure

    fig=figure('Name',fln,'NumberTitle','off');
    ax = gca(fig);
%     ax.ColorOrder=[0.7 0.7 0.7; 0.7 0.7 0.7;]; %WT r5
 ax.ColorOrder=[0.6 0.87 0.87; 0.6 0.87 0.87];%S421D r5
%     ax.ColorOrder=[0.729 0.831 0.957; 0.35 0.4 0.45; 0.549 0.667 0.745]; %WT lyso
%  ax.ColorOrder=[0.729 0.831 0.957; 0.729 0.831 0.957]; %WT lyso
%     ax.ColorOrder=[0.757 0.867 0.776;0.757 0.867 0.776];%S421D lyso
%     ax.ColorOrder=[0 0.9 0.9; 0 0.4 0.4; 0 0.75 0.75; 0 0.55 0.55];%S421D r5
%     ax.ColorOrder=[0.729 0.831 0.957; 0.35 0.4 0.45; 0.549 0.667 0.745]; %WT lyso
%     ax.ColorOrder=[0.757 0.867 0.776; 0.35 0.4 0.35; 0.522 0.74 0.549];%S421D lyso

    hold on
    cellfun(@(X,Y) plot(ax, X, Y), xk,yk) %Uncomment to plot xy positions of trajectories
%     cellfun(@(X,Y) plot(ax, X, Y), tim, pos_tim) %Uncomment to plot
%     position over time
    publication_fig(0,0,1);
%     xlabel('xk (\mum)'); ylabel('yk (\mum)');
    xlabel('Time (s)'); ylabel('Position (\mum)');
%     axis equal
    hold off
    
    
    fl_nm=[fls(k).name, nm];
%     saveas(gca,fl_nm,'png')
    
    clear res dat_pos_tim xk yk rg_2 xk_yk tim pos_tim; %Emily added tim and pos_tim Mar 28/22
    kf=0;
end
    
end

% clear ab xk yk pos kf;



toc
