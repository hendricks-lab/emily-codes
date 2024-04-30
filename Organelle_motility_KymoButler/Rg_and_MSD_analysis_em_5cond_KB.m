% Code #2 in the per trajectory analysis.
%Steps to run this code:
% 1.Update the addpaths (lines 10-15) and the directories (lines 29-91), and your DT.
% 2.Update the directory for saving figures (if you want to save them),
% lines (275-277).
% 3.Run the code!

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/breakyaxis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
% addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
% %This is not needed in this code yet.

min_lr=0.8; %BDNF=0.5, Lyso=0.8, Mito=0.6
min_disp=0;

colour_18Q=[0.5 0.5 0.5];
colour_30Q=[0.5 0.5 0.5];
% colour_30Q=[1 0 0];
colour_45Q=[0.75 0 0];
colour_65Q=[0.5 0 0];
colour_81Q=[0.25 0 0];

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good

for k_choose = 1:4
if k_choose == 1    % 30Q 
    kcol=colour_30Q;%[0.5 0.5 0.5];
    cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_proj_along_kymo_test/');
    fl='*.mat';
    an=0;
    DT=0.12; %delta t, the exposure time-modify this based on your experiment
    save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_pos/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_proj_along_kymo_test/Results/';
    kf=0;
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=1;
elseif k_choose == 2   % 45Q 
    kcol=colour_45Q;%[0.5 0.5 0.5];
%     cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_mats/');
    cd('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_proj_along_kymo_test/');
    fl='*.mat';
    an=0;
    DT=0.12; %delta t, the exposure time-modify this based on your experiment
%     save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_pos/';
    save_dir='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_pos/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_proj_along_kymo_test/Results/';
    kf=0;
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=0.2;
    
elseif k_choose == 3    % 65Q 
    kcol = colour_65Q;%[0.290000 0.330000 0.130000];
%     cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_mats/');
    cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/');
    fl='*.mat';
    an=0.030;
    DT=0.12; %delta t, the exposure time-modify this based on your experiment
%     save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_pos/';
    save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_pos/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
    kf=0;
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=0.2;
elseif k_choose == 4    % 81Q 
    kcol = colour_81Q;%[0.290000 0.330000 0.130000];
%     cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/');
    cd('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_mats/');
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_proj_along_kymo_test/');
    fl='*.mat';
    an=0.030;
    DT=0.12; %delta t, the exposure time-modify this based on your experiment
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_pos/';
    save_dir='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_pos/';
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_proj_along_kymo_test/Results/';
    kf=0;
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=0.2;
    
end

fls=dir(fl);
kpr=0;
    
for k=1:numel(fls)
    load(fls(k).name);
    display(fls(k).name);
    for j=1:numel(kbpos)
%     for j=1:numel(ab)
%         disp(j)
    if (numel(kbpos(j).position)>15) %&& (mean(ab(j).medint>3000))
    kf=kf+1;
    pos{kf}=kbpos(j).position;
    tim{kf}=[1:numel(kbpos(j).position)]'.*DT;%ab(j).tk;
    position_smooth=smooth(kbpos(j).position,span,'sgolay',pwr);

    %% Radius of Gyration:
    rg_1=func_Rg_1D(kbpos(j).position);
    rg=[rg,rg_1];

%     if (numel(ab(j).position)>15) %&& (mean(ab(j).medint>3000))
%     kf=kf+1;
%     pos{kf}=ab(j).position;
%     tim{kf}=[1:numel(ab(j).position)]'.*DT;%ab(j).tk;
%     position_smooth=smooth(ab(j).position,span,'sgolay',pwr);
% 
%     %% Radius of Gyration:
%     rg_1=func_Rg_1D(ab(j).position);
%     rg=[rg,rg_1];
    
    %% Mean Squared Displacement
    [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_1d (position_smooth, DT, k_choose);
%     res= msd_fun_v1_emmessedwithit(40,position_smooth,DT,numel(position_smooth),0);
%     alpha_1=res.slp;
%     msdpts=res.msd;
%     deltat=res.timek;

    if alpha_1>0
    alpha=[alpha,alpha_1];
    delT{kf}=deltat;
    MSD{kf}=msdpts;
    dat_pos_tim{kf}=position_smooth;
    Diff_coeff=[Diff_coeff,DiffCoef];
    else
%         disp('Alpha less than or equal to 0');
    end
    
    else
    end
    clear res
    end
end
    kf
    
    dat_pos_tim=dat_pos_tim(~cellfun('isempty',dat_pos_tim));
    delT=delT(~cellfun('isempty',delT));
    MSD=MSD(~cellfun('isempty',MSD));
   
    rg_kch{k_choose}=rg;
    alpha_kch{k_choose}=alpha;
    Diff_coeff_kch{k_choose}=Diff_coeff;
    position_dat{k_choose}=dat_pos_tim;
   
    display(mean(alpha))
    
    clear delT MSD rg_hist pos tim dat_pos_tim frac_plus_time frac_minus_time x_hist_plus x_hist_minus rg_hist rg_x mean_intensity;

    position=position_dat{k_choose};
    Diffusion_coefficient=Diff_coeff_kch{k_choose};

    save([save_dir, 'position'],'position','Diffusion_coefficient');
    
    clear position Diffusion_coefficient
    
    proc_runs=[];
    proc_times=[];
    nvproc=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
end

%% Plot alpha distribution:
alpha_30Q=alpha_kch{1};
alpha_45Q=alpha_kch{2};
alpha_65Q=alpha_kch{3};
alpha_81Q=alpha_kch{4};


alpha_x=0.2:0.04:2;
alpha_hist_30Q=hist(alpha_30Q,alpha_x);
alpha_hist_45Q=hist(alpha_45Q,alpha_x);
alpha_hist_65Q=hist(alpha_65Q,alpha_x);
alpha_hist_81Q=hist(alpha_81Q,alpha_x);

figure('Name','Alpha_Distribution','NumberTitle','off'), hold on, 
bar(alpha_x, alpha_hist_30Q./sum(alpha_hist_30Q), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',0.3);
hold on,
bar(alpha_x, alpha_hist_45Q./sum(alpha_hist_45Q), 'BarWidth',1,'FaceColor',[0, 0.5, 0.75],'edgecolor',[0, 0.5, 0.75],'Facealpha',0.3);
hold on,
bar(alpha_x, alpha_hist_65Q./sum(alpha_hist_65Q), 'BarWidth',1,'FaceColor',[0, 0.3, 0.75],'edgecolor',[0, 0.3, 0.75],'Facealpha',0.3);
hold on,
bar(alpha_x, alpha_hist_81Q./sum(alpha_hist_81Q), 'BarWidth',1,'FaceColor',[0, 0.1, 0.75],'edgecolor',[0, 0.1, 0.75],'Facealpha',0.3);
hold on,
xlabel('\alpha'); ylabel('Number of Trajectories');
xlim([0 2.1]);
publication_fig(0,0,1);

%% Plot Rg distribution:
% rg_18Q=rg_kch{1};
rg_30Q=rg_kch{1};
rg_45Q=rg_kch{2};
rg_65Q=rg_kch{3};
rg_81Q=rg_kch{4};

%rg_x=0.01:0.03:2.5;
rg_x=0:0.03:2.5;
% rg_hist_18Q=hist(rg_18Q,rg_x);
rg_hist_30Q=hist(rg_30Q,rg_x);
rg_hist_45Q=hist(rg_45Q,rg_x);
rg_hist_65Q=hist(rg_65Q,rg_x);
rg_hist_81Q=hist(rg_81Q,rg_x);
figure('Name','Rg_Distribution','NumberTitle','off'), hold on, hold on, 
bar(rg_x, rg_hist_30Q./sum(rg_hist_30Q), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',0.2);
hold on,
bar(rg_x, rg_hist_45Q./sum(rg_hist_45Q), 'BarWidth',1,'FaceColor',[0, 0.5, 0.75],'edgecolor',[0, 0.5, 0.75],'Facealpha',0.2);
hold on,
bar(rg_x, rg_hist_65Q./sum(rg_hist_65Q), 'BarWidth',1,'FaceColor',[0, 0.3, 0.75],'edgecolor',[0, 0.3, 0.75],'Facealpha',0.2);
hold on,
bar(rg_x, rg_hist_81Q./sum(rg_hist_81Q), 'BarWidth',1,'FaceColor',[0, 0.1, 0.75],'edgecolor',[0, 0.1, 0.75],'Facealpha',0.2);
hold on,
xlabel('Radius of Gyration (\mum)'); ylabel('Number of Trajectories');
xlim([0 2.1]);
publication_fig(0,0,1);


figure('Name','Rg_per_trajectory_boxplot','NumberTitle','off'), hold on,
h2=notBoxPlot(rg_kch{1},1,'style','line');
set(h2.data,'MarkerFaceColor',colour_30Q,'MarkerEdgeColor',colour_30Q);
h3=notBoxPlot(rg_kch{2},2,'style','line');
set(h3.data,'MarkerFaceColor',colour_45Q,'MarkerEdgeColor',colour_45Q);
h4=notBoxPlot(rg_kch{3},3,'style','line');
set(h4.data,'MarkerFaceColor',colour_65Q,'MarkerEdgeColor',colour_65Q);
h5=notBoxPlot(rg_kch{4},4,'style','line');
set(h5.data,'MarkerFaceColor',colour_81Q,'MarkerEdgeColor',colour_81Q);
set(gca,'Ygrid','on');
set(gca,'FontName','Arial','FontSize',20,'LineWidth',2,'XColor',[0 0 0],...
    'XTick',[1 2 3 4 5],'XTickLabel',{'18Q','30Q','45Q','65Q','81Q'},...
    'YColor',[0 0 0]);
ylabel('Radius of Gyration (\mum)');
xlim([0 6]);
ylim([-0.5 9]);
publication_fig(0,0,1)


figure('Name','Rg_per_trajectory_histogram','NumberTitle','off'), hold on,
histogram(rg_kch{1}, 'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.2), hold on,
histogram(rg_kch{2}, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.2), hold on,
histogram(rg_kch{3}, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
histogram(rg_kch{4}, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);


% fileprefix='20230807_isoHD_bdnf_'
% 
% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/lyso/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end


