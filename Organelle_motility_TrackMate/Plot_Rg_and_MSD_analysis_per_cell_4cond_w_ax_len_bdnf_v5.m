%% Analyze and Plot Rg and MSD per cell
%% Abdullah R. Chaudhary (modified by Emily N. P. Prowse)

%Code #3 to run (after analyze trackmate motility and analyze dispersion).

%Steps for running this code:
% 1.Update the addpaths (lines 17-21) and directories (lines 40-129) and variable DT.
% 2.Update the directory for saving figures (if needed) (lines 376-378)
% 3.If you have less than 5 conditions, update the number of k_choose,
% comment out the sections for the unused conditions, and update the violin
% plot axis labels and colours manually.
% 4. Run the code!


close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
% %None of this directory is used in the current version of the code,but if
% statistical testing parts are updated they will be. 

plus_processive_runs=[];

%Setting the colour code for each condition (because it's so annoying to do
%in each plot)
% polyQ length
cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_45Q=cmap(20,:);%[0.75 0 0];
colour_65Q=cmap(25,:);%[0.5 0 0];
colour_81Q=cmap(30,:);%[0.25 0 0];
htt_30Q="30Q";
htt_45Q="45Q";
htt_65Q="65Q";
htt_81Q="81Q";

%stress condition
% cmap = colormap(hot(100));
% colour_30Q=[0.5 0.5 0.5];
% colour_65Q=cmap(30,:);%[0.5 0 0];
% cmap2= colormap(spring(100));
% colour_45Q=[1.0000 0.4444 0.5556]%cmap2(50,:);%[0.75 0 0];
% colour_81Q=cmap2(30,:);%[0.25 0 0];
% 
% htt_30Q="30Q";
% htt_45Q="30Q+IFg";
% htt_65Q="81Q";
% htt_81Q="81Q+IFg";

fileprefix='20231012_isoHD_mito_proj_';

for k_choose = 1:4

if k_choose == 1    % 30Q
    col1=colour_30Q;
    %Lyso projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_proj_along_kymo_test/Results/';
    %Lyso
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_dir/';
%30Q kymobutler lyso
%     pth1='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_dir/';
%     %bdnf
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_dir/';
    %30Q kymobutler bdnf
%     pth1='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_dir/';
%BDNF projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_proj_along_kymo_test/Results/';
    %Mito
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_dir/';
    %30Q kymobutler mito
%     pth1='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_dir/';
%Mito projected along kymograph
    pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_proj_along_kymo_test/Results/';
    pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_proj_along_kymo_test/Results/';
    pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_proj_along_kymo_test/Results/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
%     fl10='axon_length_per_cell.mat';
%     fl11='alphas.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];

elseif k_choose == 2   % 45Q
    col1=colour_45Q;
    %45Q Lyso projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_proj_along_kymo_test/Results/';
    %Lyso IFg projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_proj_along_kymo_test/Results/';
    %Lyso
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_dir/';
    %45Q Lyso kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_dir/';
    %kymobutler bdnf+IFg
%     pth1='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_dir/';
    %Bdnf
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_dir/';
% %BDNF projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_proj_along_kymo_test/Results/';
% %BDNF IFg projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_proj_along_kymo_test/Results/';
 %Mito
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_dir/';
    %45Q kymobutler mito
%     pth1='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_dir/';
%45Q Mito projected along kymograph
    pth1='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_proj_along_kymo_test/Results/';
    pth2='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_proj_along_kymo_test/Results/';
    pth3='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_proj_along_kymo_test/Results/';
    %30Q bdnf+IFg
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_dir/';
    %45Q BDNF kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_dir/';
    %kymobutler bdnf+IFg
%     pth1='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_dir/';
    %30Q lyso+IFg
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_dir/';
   
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
%     fl10='axon_length_per_cell.mat';
%     fl11='alphas.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=1;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
    
elseif k_choose == 3   % 65Q
    col1=colour_65Q;
    % 65Q Lyso projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_proj_along_kymo_test/Results/';
    %81Q Lyso projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
    % 65Q Lyso
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_dir/';
    %65Q Lyso kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_dir/';
    %81Q Lyso kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_dir/';
    % 65Q Bdnf
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_dir/';
    % 65Q bdnf projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_proj_along_kymo_test/Results/';
    %81Q bdnf projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
    %Mito
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_dir/';
    %65Q mito kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_dir/';
    % 65Q Mito projected along kymograph
    pth1='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_proj_along_kymo_test/Results/';
    pth2='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_proj_along_kymo_test/Results/';
    pth3='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_proj_along_kymo_test/Results/';
    %81Q bdnf
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_dir/';
    %65Q BDNF kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_dir/';
    %81Q BDNF kymobutler
%     pth1='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_dir/';
     %81Q lyso
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
%     fl10='axon_length_per_cell.mat';
%     fl11='alphas.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=1.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
    
 elseif k_choose == 4   % 81Q
    col1=colour_81Q;
    %Lyso projected along kymograph    
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_proj_along_kymo_test/Results/';
%     %Lyso IFg projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_proj_along_kymo_test/Results/';
    % 81Q bdnf projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_proj_along_kymo_test/Results/';
    %81Q bdnf IFg projected along kymograph
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_proj_along_kymo_test/Results/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_proj_along_kymo_test/Results/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_proj_along_kymo_test/Results/';
%Mito projected along kymograph    
    pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_proj_along_kymo_test/Results/';
    pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_proj_along_kymo_test/Results/';
    pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_proj_along_kymo_test/Results/';
    %Lyso
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_dir/';
%   Kymobutler lyso 81Q
%     pth1='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_dir/';
    %kymobutler 81Q lyso+IFg
%     pth1='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_dir/';
%     %Bdnf
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_dir/';
 %Mito
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_dir/';
%   Kymobutler Mito 81Q
%     pth1='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_dir/';
    %81Q bdnf+IFg
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_dir/';
%   Kymobutler BDNF 81Q
%     pth1='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_dir/';
    %kymobutler 81Q bdnf+IFg
%     pth1='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_rg/';
%     pth2='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_msd/';
%     pth3='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_dir/';
    %81Q lyso+IFg
%     pth1='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_rg/';
%     pth2='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_msd/';
%     pth3='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_dir/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
%     fl10='axon_length_per_cell.mat';
%     fl11='alphas.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=1.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
end

load(fullfile(pth1,fl1));
load(fullfile(pth2,fl2));
load(fullfile(pth3,fl3)); %Emily added
load(fullfile(pth3,fl4)); %Emily added 20220927
load(fullfile(pth3,fl5)); %Emily added 20220927
load(fullfile(pth3,fl6)); %Emily added 20220927
load(fullfile(pth3,fl7)); %Emily added 20221020
load(fullfile(pth3,fl8)); %Emily added 20221020
load(fullfile(pth3,fl9));
% load(fullfile(pth3, fl10)); %Emily added 20230208
% load(fullfile(pth2, fl11)); %Emily added 20230208

%Emily adding in a longer minimum run length version 20221024
proc_frac_t{k_choose}=t_proc_per_cell;
plus_proc_frac_t{k_choose}=t_proc_out_per_cell;
minus_proc_frac_t{k_choose}=t_proc_in_per_cell;
diff_frac_t{k_choose}=t_diff_per_cell;
plus_diff_frac_t{k_choose}=t_diff_out_per_cell; %Emily added 20221020
minus_diff_frac_t{k_choose}=t_diff_in_per_cell; %Emily added 20221020
stat_frac_t{k_choose}=t_stat_per_cell;
% axon_length{k_choose}=axon_length_per_cell; %Emily added 20230208
% alpha_2{k_choose}=cell2mat(alpha_2d);
[timek2,msd2] = average_MSD_per_cell_emmodified4kchoose(res1,col1,k_choose);

% MSD Analysis:
for ik=1:numel(res1)
    timek=res1{ik}.timek;
    msd=res1{ik}.msd;
    logt=res1{ik}.logtk;
    logmsd=res1{ik}.log_msd;
    alph=res1{ik}.slp;
    if isnan(msd) %Emily added because the msd had some cells without any data
        k_choose %1 for 18Q 2 for 30Q, helps identify the empty cell
        ik
        warning='Empty cell in msd analysis'
    else

    % Figure(1)

    figure(k_choose*1e1), hold on, 
    p1=plot(timek,msd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Time interval (sec)'); ylabel('MSD (\mum^2)'); 
    publication_fig(0,0,1);
    
    % Figure(1)
    figure(k_choose*1e3), hold on, 
    p1=plot(logt,logmsd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Log[Time interval (sec)]'); ylabel('Log[MSD (\mum^2)]'); 
    publication_fig(0,0,1);
    
    msd_all=[msd_all,msd];
    slp_all=[slp_all,alph];
    
    rg_per_cell=rg_all{ik};
    if isempty(rg_per_cell) %Emily added because the Rg had some cells without any data
        k_choose %1 for 18Q 2 for 30Q, helps identify the empty cell
        ik
        warning='Empty cell in Rg analysis'
        %clear rg_per_cell
    else
    [mean_rg,pci]=mle(rg_per_cell,'distribution','exp'); 
    mean_rg_all=[mean_rg_all,mean_rg];
    
    %Emily plotting the max likelihood estimation of the rg and to figure out what's wrong
%     rg_x=0:0.03:2.5;
%     rg_hist_per_cell=hist(rg_per_cell,rg_x);
%     rg_hist_mean=hist(mean_rg_all,rg_x);
    
%     figure(k_choose*8), hold on, 
%     bar(rg_x, rg_hist_mean./sum(rg_hist_mean), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',1);
%     hold on,
%     stairs(rg_x, rg_hist_per_cell./sum(rg_hist_per_cell), 'LineWidth',2,'Color','k');
%     hold on, 
%     xlabel('Radius of Gyration (\mum)'); ylabel('Number of Trajectories');
%     xlim([0 2.1]);
%     publication_fig(0,0,1);

%     figure(k_choose*9), hold on,
%     scatter(1:numel(rg_per_cell),rg_per_cell,'r')
%     plot(1:numel(mean_rg_all),mean_rg_all,'g');
    clear rg_per_cell
    end
    end
%     [mean_rg,pci]=mle(rg_per_cell,'distribution','exp'); %original
%     %Abdullah's code but for some reason my data had some NaNs in it so
%     %we had to add the if: else statement. 
%     mean_rg_all=[mean_rg_all,mean_rg];
%     clear rg_per_cell
end

msd_data{k_choose}=msd_all';
rg_data{k_choose}=mean_rg_all';
timek_f{k_choose}=timek2;
MSD_f{k_choose}=msd2;
alph_f{k_choose}=slp_all; %slopes from the fits to the MSD curve

clear msd_all mean_rg_all rg_per_cell timek2 msd2 slp_all outward_runs_per_cell;
%rmpath(pth1)
%rmpath(pth2)
end
%This is required when there are some NaNs in the alpha/rg code
% axon_length{1}=axon_length{1}(~isnan(alph_f{1}));
% axon_length{2}=axon_length{2}(~isnan(alph_f{2}));
% axon_length{3}=axon_length{3}(~isnan(alph_f{3}));
% axon_length{4}=axon_length{4}(~isnan(alph_f{4}));

plus_proc_frac_t{1}=plus_proc_frac_t{1}(~isnan(alph_f{1}));
plus_proc_frac_t{2}=plus_proc_frac_t{2}(~isnan(alph_f{2}));
plus_proc_frac_t{3}=plus_proc_frac_t{3}(~isnan(alph_f{3}));
plus_proc_frac_t{4}=plus_proc_frac_t{4}(~isnan(alph_f{4}));

%Per cell analysis
figure('Name','Alpha_per_cell_violin_and_box','NumberTitle','off'), hold on,
violin(alph_f, 'xlabel',{'30Q','45Q','65Q','81Q'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% violin(alph_f, 'xlabel',{'30Q','30Q+IFg','81Q','81Q+IFg'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(alph_f{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerSize', 5);
h2=notBoxPlot_nodotorline(alph_f{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(alph_f{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerSize', 5);
h4=notBoxPlot_nodotorline(alph_f{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerSize', 5);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('\alpha');
% ylim([0.5 2.25]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
xticks([1 2 3 4]);
xticklabels({'30Q', '45Q','65Q','81Q'});
% xticklabels({'30Q', '30Q+IFγ','81Q','81Q+IFγ'});
box("off");
set(gca,'FontSize',40);

a_fraction_change{1}=(mean(alph_f{1},'omitnan')-mean(alph_f{1},'omitnan'))/mean(alph_f{1},'omitnan');
a_fraction_change{2}=(mean(alph_f{2},'omitnan')-mean(alph_f{1},'omitnan'))/mean(alph_f{1},'omitnan');
a_fraction_change{3}=(mean(alph_f{3},'omitnan')-mean(alph_f{1},'omitnan'))/mean(alph_f{1},'omitnan');
a_fraction_change{4}=(mean(alph_f{4},'omitnan')-mean(alph_f{1},'omitnan'))/mean(alph_f{1},'omitnan');

a_fold_change{1}=mean(alph_f{1},'omitnan')/mean(alph_f{1},'omitnan');
a_fold_change{2}=mean(alph_f{2},'omitnan')/mean(alph_f{1},'omitnan');
a_fold_change{3}=mean(alph_f{3},'omitnan')/mean(alph_f{1},'omitnan');
a_fold_change{4}=mean(alph_f{4},'omitnan')/mean(alph_f{1},'omitnan');

figure('Name','Alpha_per_cell_means_fraction_change','NumberTitle','off'), hold on,
plot(0.5,a_fraction_change{1},'_','MarkerSize',80,'MarkerEdgeColor', colour_30Q, 'LineWidth', 10);
plot(1.5,a_fraction_change{2},'_','MarkerSize',80,'MarkerEdgeColor', colour_45Q, 'LineWidth', 10);
plot(2.5,a_fraction_change{3},'_','MarkerSize',80,'MarkerEdgeColor', colour_65Q, 'LineWidth', 10);
plot(3.5,a_fraction_change{4},'_','MarkerSize',80,'MarkerEdgeColor', colour_81Q, 'LineWidth', 10);
plot([0 4],[0 0],'k-','LineWidth',0.5);
set(gca, "FontSize", 30);
xlim([0 4]);
ylabel('Δ');
h = gca;
h.XAxis.Visible = 'off';
pbaspect([4 1 1]);

figure('Name','Rg_per_cell_violin_and_box','NumberTitle','off'), hold on,
violin(rg_data, 'xlabel',{'30Q','45Q','65Q','81Q'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% violin(rg_data, 'xlabel',{'30Q','30Q+IFg','81Q','81Q+IFg'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(rg_data{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(rg_data{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(rg_data{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerSize', 5);
h4=notBoxPlot_nodotorline(rg_data{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerSize', 5);
% xlim([0.5 5.5]);
% ylim([-0.25 2]);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Radius of Gyration (\mum)');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')
xticks([1 2 3 4]);
xticklabels({'30Q', '45Q','65Q','81Q'});
% xticklabels({'30Q', '30Q+IFγ','81Q','81Q+IFγ'});
box("off");
set(gca,'FontSize',40);

rg_fraction_change{1}=(mean(rg_data{1},'omitnan')-mean(rg_data{1},'omitnan'))/mean(rg_data{1},'omitnan');
rg_fraction_change{2}=(mean(rg_data{2},'omitnan')-mean(rg_data{1},'omitnan'))/mean(rg_data{1},'omitnan');
rg_fraction_change{3}=(mean(rg_data{3},'omitnan')-mean(rg_data{1},'omitnan'))/mean(rg_data{1},'omitnan');
rg_fraction_change{4}=(mean(rg_data{4},'omitnan')-mean(rg_data{1},'omitnan'))/mean(rg_data{1},'omitnan');

rg_fold_change{1}=mean(rg_data{1},'omitnan')/mean(rg_data{1},'omitnan');
rg_fold_change{2}=mean(rg_data{2},'omitnan')/mean(rg_data{1},'omitnan');
rg_fold_change{3}=mean(rg_data{3},'omitnan')/mean(rg_data{1},'omitnan');
rg_fold_change{4}=mean(rg_data{4},'omitnan')/mean(rg_data{1},'omitnan');

figure('Name','Rg_per_cell_means_fraction_change','NumberTitle','off'), hold on,
plot(0.5,rg_fraction_change{1},'_','MarkerSize',80,'MarkerEdgeColor', colour_30Q, 'LineWidth', 10);
plot(1.5,rg_fraction_change{2},'_','MarkerSize',80,'MarkerEdgeColor', colour_45Q, 'LineWidth', 10);
plot(2.5,rg_fraction_change{3},'_','MarkerSize',80,'MarkerEdgeColor', colour_65Q, 'LineWidth', 10);
plot(3.5,rg_fraction_change{4},'_','MarkerSize',80,'MarkerEdgeColor', colour_81Q, 'LineWidth', 10);
plot([0 4],[0 0],'k-','LineWidth',0.5);
set(gca, "FontSize", 30);
xlim([0 4]);
ylabel('Δ');
h = gca;
h.XAxis.Visible = 'off';
pbaspect([4 1 1]);

figure('Name','Directional bias run length box and violin','NumberTitle','off'), hold on,
violin(plus_proc_frac_t, 'xlabel',{'30Q','45Q','65Q','81Q'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% violin(plus_proc_frac_t, 'xlabel',{'30Q','30Q+IFg','81Q','81Q+IFg'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(plus_proc_frac_t{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(plus_proc_frac_t{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(plus_proc_frac_t{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerSize', 5);
h4=notBoxPlot_nodotorline(plus_proc_frac_t{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerSize', 5);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Fraction of time of processive runs moving outward');
% ylim([0.2 0.75]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

% %Emily added for testing if min_lr was set correctly 20221020
% figure('Name','Diff_runs_in_out_box','NumberTitle','off'), hold on,
% violin(plus_diff_frac_t, 'xlabel',{'30Q','45Q','65Q','81Q'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% % violin(plus_diff_frac_t, 'xlabel',{'30Q','30Q+IFg','81Q','81Q+IFg'},'facecolor',[colour_30Q;colour_45Q;colour_65Q;colour_81Q],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% h5=notBoxPlot_nodotorline(plus_diff_frac_t{1},1,'style','line');
% set(h5.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q, 'MarkerSize', 5);
% h6=notBoxPlot_nodotorline(plus_diff_frac_t{2},2,'style','line');
% set(h6.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerSize', 5);
% h7=notBoxPlot_nodotorline(plus_diff_frac_t{3},3,'style','line');
% set(h7.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerSize', 5);
% h8=notBoxPlot_nodotorline(plus_diff_frac_t{4},4,'style','line');
% set(h8.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerSize', 5);
% ylabel('Fraction of time of diffusive runs moving outwards');
% publication_fig(0,0,1);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off');

%Removing NaNs in motility analysis for BDNF condition
plus_proc_ft{1}=plus_proc_frac_t{1}(~isnan(plus_proc_frac_t{1}));
plus_proc_ft{2}=plus_proc_frac_t{2}(~isnan(plus_proc_frac_t{2}));
plus_proc_ft{3}=plus_proc_frac_t{3}(~isnan(plus_proc_frac_t{3}));
plus_proc_ft{4}=plus_proc_frac_t{4}(~isnan(plus_proc_frac_t{4}));
minus_proc_ft{1}=minus_proc_frac_t{1}(~isnan(minus_proc_frac_t{1}));
minus_proc_ft{2}=minus_proc_frac_t{2}(~isnan(minus_proc_frac_t{2}));
minus_proc_ft{3}=minus_proc_frac_t{3}(~isnan(minus_proc_frac_t{3}));
minus_proc_ft{4}=minus_proc_frac_t{4}(~isnan(minus_proc_frac_t{4}));
proc_ft{1}=proc_frac_t{1}(~isnan(proc_frac_t{1}));
proc_ft{2}=proc_frac_t{2}(~isnan(proc_frac_t{2}));
proc_ft{3}=proc_frac_t{3}(~isnan(proc_frac_t{3}));
proc_ft{4}=proc_frac_t{4}(~isnan(proc_frac_t{4}));
diff_ft{1}=diff_frac_t{1}(~isnan(diff_frac_t{1}));
diff_ft{2}=diff_frac_t{2}(~isnan(diff_frac_t{2}));
diff_ft{3}=diff_frac_t{3}(~isnan(diff_frac_t{3}));
diff_ft{4}=diff_frac_t{4}(~isnan(diff_frac_t{4}));
stat_ft{1}=stat_frac_t{1}(~isnan(stat_frac_t{1}));
stat_ft{2}=stat_frac_t{2}(~isnan(stat_frac_t{2}));
stat_ft{3}=stat_frac_t{3}(~isnan(stat_frac_t{3}));
stat_ft{4}=stat_frac_t{4}(~isnan(stat_frac_t{4}));

mot_frac{1}=proc_ft{1};
mot_frac{2}=proc_ft{2};
mot_frac{3}=proc_ft{3};
mot_frac{4}=proc_ft{4};
mot_frac{5}=diff_ft{1};
mot_frac{6}=diff_ft{2};
mot_frac{7}=diff_ft{3};
mot_frac{8}=diff_ft{4};
mot_frac{9}=stat_ft{1};
mot_frac{10}=stat_ft{2};
mot_frac{11}=stat_ft{3};
mot_frac{12}=stat_ft{4};

%Mean all conditions
plus_runs_30Q=mean(plus_proc_ft{1});
plus_runs_45Q=mean(plus_proc_ft{2});
plus_runs_65Q=mean(plus_proc_ft{3});
plus_runs_81Q=mean(plus_proc_ft{4});
minus_runs_30Q=mean(minus_proc_ft{1});
minus_runs_45Q=mean(minus_proc_ft{2});
minus_runs_65Q=mean(minus_proc_ft{3});
minus_runs_81Q=mean(minus_proc_ft{4});
proc_runs_30Q=mean(proc_ft{1});
proc_runs_45Q=mean(proc_ft{2});
proc_runs_65Q=mean(proc_ft{3});
proc_runs_81Q=mean(proc_ft{4});
diff_frac_30Q=mean(diff_ft{1});
diff_frac_45Q=mean(diff_ft{2});
diff_frac_65Q=mean(diff_ft{3});
diff_frac_81Q=mean(diff_ft{4});
stat_frac_30Q=mean(stat_ft{1});
stat_frac_45Q=mean(stat_ft{2});
stat_frac_65Q=mean(stat_ft{3});
stat_frac_81Q=mean(stat_ft{4});

%SEM for directional bias parameters
error_proc_30Q=std(proc_ft{1})/(sqrt(length(proc_ft{1})));
error_diff_30Q=std(diff_ft{1})/(sqrt(length(diff_ft{1})));
error_stat_30Q=std(stat_ft{1})/(sqrt(length(stat_ft{1})));
error_proc_45Q=std(proc_ft{2})/(sqrt(length(proc_ft{2})));
error_diff_45Q=std(diff_ft{2})/(sqrt(length(diff_ft{2})));
error_stat_45Q=std(stat_ft{2})/(sqrt(length(stat_ft{2})));
error_proc_65Q=std(proc_ft{3})/(sqrt(length(proc_ft{3})));
error_diff_65Q=std(diff_ft{3})/(sqrt(length(diff_ft{3})));
error_stat_65Q=std(stat_ft{3})/(sqrt(length(stat_ft{3})));
error_proc_81Q=std(proc_ft{4})/(sqrt(length(proc_ft{4})));
error_diff_81Q=std(diff_ft{4})/(sqrt(length(diff_ft{4})));
error_stat_81Q=std(stat_ft{4})/(sqrt(length(stat_ft{4})));

%creating a variable with all the errors for all the conditions
all_err_bars=[ error_proc_30Q error_proc_45Q error_proc_65Q error_proc_81Q; error_diff_30Q error_diff_45Q error_diff_65Q error_diff_81Q; error_stat_30Q error_stat_45Q error_stat_65Q error_stat_81Q];
% 
values=[ proc_runs_30Q proc_runs_45Q proc_runs_65Q proc_runs_81Q; diff_frac_30Q diff_frac_45Q diff_frac_65Q diff_frac_81Q; stat_frac_30Q stat_frac_45Q stat_frac_65Q stat_frac_81Q];
% all_err_bars=[error_proc_18Q error_proc_30Q error_proc_45Q error_proc_81Q; error_diff_18Q error_diff_30Q error_diff_45Q error_diff_81Q; error_stat_18Q error_stat_30Q error_stat_45Q error_stat_81Q];

% values=[proc_runs_18Q proc_runs_30Q proc_runs_45Q proc_runs_81Q; diff_frac_18Q diff_frac_30Q diff_frac_45Q diff_frac_81Q;stat_frac_18Q stat_frac_30Q stat_frac_45Q stat_frac_81Q];
figure('Name','Bar plot proc stat diff','NumberTitle','off'), hold on,
b2=bar(values,'grouped');
b2(1).FaceColor=colour_30Q;
b2(2).FaceColor=colour_45Q;
b2(3).FaceColor=colour_65Q;
b2(4).FaceColor=colour_81Q;
b2(1).FaceAlpha=0.5;
b2(2).FaceAlpha=0.5;
b2(3).FaceAlpha=0.5;
b2(4).FaceAlpha=0.5;
hold on
scatter(0.7273-0.02+0.05*rand(numel(mot_frac{1}),1),mot_frac{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(0.9091-0.02+0.05*rand(numel(mot_frac{2}),1),mot_frac{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.0909-0.02+0.05*rand(numel(mot_frac{3}),1),mot_frac{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.2727-0.02+0.05*rand(numel(mot_frac{4}),1),mot_frac{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.7273-0.02+0.05*rand(numel(mot_frac{5}),1),mot_frac{5},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1.9091-0.02+0.05*rand(numel(mot_frac{6}),1),mot_frac{6},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.0909-0.02+0.05*rand(numel(mot_frac{7}),1),mot_frac{7},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.2727-0.02+0.05*rand(numel(mot_frac{8}),1),mot_frac{8},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.7273-0.02+0.05*rand(numel(mot_frac{9}),1),mot_frac{9},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2.9091-0.02+0.05*rand(numel(mot_frac{10}),1),mot_frac{10},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3.0909-0.02+0.05*rand(numel(mot_frac{11}),1),mot_frac{11},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3.2727-0.02+0.05*rand(numel(mot_frac{12}),1),mot_frac{12},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
[ngroups,nbars] = size(values);
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, values(:,i), all_err_bars(:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
ylabel('Fraction of cargoes');
xticks([1 2 3]);
xticklabels({'Processive', 'Diffusive','Stationary'});
ylabel('Fraction of runs');
publication_fig(0,0,1)
box("off");
set(gca,'FontSize',40);


%Statistical testing section

%Rg t test
[rg_30Q_45Q,p_30Q_45Q_rg]=ttest2(rg_data{1}, rg_data{2});
[rg_30Q_65Q,p_30Q_65Q_rg]=ttest2(rg_data{1}, rg_data{3});
[rg_30Q_81Q,p_30Q_81Q_rg]=ttest2(rg_data{1}, rg_data{4});

%Rg one-way ANOVA
rg_values=[rg_data{1};rg_data{2};rg_data{3};rg_data{4}];
rg_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(rg_data{1}) numel(rg_data{2}) numel(rg_data{3}) numel(rg_data{4})])]';
rg_data_tbl=table(rg_values,rg_groups);

[rg_p,rg_tbl_anova,rg_stats]=anova1(rg_values,rg_groups,"on");

[rg_results,~,~,gnames]=multcompare(rg_stats);

rg_tbl = array2table(rg_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
rg_tbl.("Group A") = gnames(rg_tbl.("Group A"));
rg_tbl.("Group B") = gnames(rg_tbl.("Group B"));

% %Alpha t test
[alpha_30Q_45Q,p_30Q_45Q_alpha]=ttest2(alph_f{1}, alph_f{2});
[alpha_30Q_65Q,p_30Q_65Q_alpha]=ttest2(alph_f{1}, alph_f{3});
[alpha_30Q_81Q,p_30Q_81Q_alpha]=ttest2(alph_f{1}, alph_f{4});

alph_values=[alph_f{1}';alph_f{2}';alph_f{3}';alph_f{4}'];
alph_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(alph_f{1}) numel(alph_f{2}) numel(alph_f{3}) numel(alph_f{4})])]';
alph_data=table(alph_values,alph_groups);

[alph_p,alph_tbl_anova,alph_stats]=anova1(alph_values,alph_groups,"on");

[alph_results,~,~,gnames]=multcompare(alph_stats);

alph_tbl = array2table(alph_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
alph_tbl.("Group A") = gnames(alph_tbl.("Group A"));
alph_tbl.("Group B") = gnames(alph_tbl.("Group B"));



% %Directional Bias t test

[dirbias_30Q_45Q,p_30Q_45Q_dirbias]=ttest2(plus_proc_ft{1}, plus_proc_ft{2});
[dirbias_30Q_65Q,p_30Q_65Q_dirbias]=ttest2(plus_proc_ft{1}, plus_proc_ft{3});
[dirbias_30Q_81Q,p_30Q_81Q_dirbias]=ttest2(plus_proc_ft{1}, plus_proc_ft{4});

dirbias_values=[plus_proc_ft{1};plus_proc_ft{2};plus_proc_ft{3};plus_proc_ft{4}];
dirbias_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(plus_proc_ft{1}) numel(plus_proc_ft{2}) numel(plus_proc_ft{3}) numel(plus_proc_ft{4})])]';
dirbias_data=table(dirbias_values,dirbias_groups);

[dirbias_p,dirbias_tbl_anova,dirbias_stats]=anova1(dirbias_values,dirbias_groups,"on");

[dirbias_results,m,h,gnames]=multcompare(dirbias_stats);

dirbias_tbl = array2table(dirbias_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
dirbias_tbl.("Group A") = gnames(dirbias_tbl.("Group A"));
dirbias_tbl.("Group B") = gnames(dirbias_tbl.("Group B"));

%Fraction proc/diff/stat stats

[t_30Q_45Q_proc,p_30Q_45Q_proc]=ttest2(proc_ft{1}, proc_ft{2});
[t_30Q_65Q_proc,p_30Q_65Q_proc]=ttest2(proc_ft{1}, proc_ft{3});
[t_30Q_81Q_proc,p_30Q_81Q_proc]=ttest2(proc_ft{1}, proc_ft{4});

[t_30Q_45Q_diff,p_30Q_45Q_diff]=ttest2(diff_ft{1}, diff_ft{2});
[t_30Q_65Q_diff,p_30Q_65Q_diff]=ttest2(diff_ft{1}, diff_ft{3});
[t_30Q_81Q_diff,p_30Q_81Q_diff]=ttest2(diff_ft{1}, diff_ft{4});

[t_30Q_45Q_stat,p_30Q_45Q_stat]=ttest2(stat_ft{1}, stat_ft{2});
[t_30Q_65Q_stat,p_30Q_65Q_stat]=ttest2(stat_ft{1}, stat_ft{3});
[t_30Q_81Q_stat,p_30Q_81Q_stat]=ttest2(stat_ft{1}, stat_ft{4});

proc="proc";
diff="diff";
stat="stat";

proc_diff_stat_all=[proc_ft{1}; diff_ft{1}; stat_ft{1}; proc_ft{2}; diff_ft{2}; stat_ft{2}; proc_ft{3}; diff_ft{3}; stat_ft{3}; proc_ft{4}; diff_ft{4}; stat_ft{4}];
cell_type=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(proc_ft{1})+numel(diff_ft{1})+numel(stat_ft{1}) numel(proc_ft{2})+numel(diff_ft{2})+numel(stat_ft{2}) numel(proc_ft{3})+numel(diff_ft{3})+numel(stat_ft{3}) numel(proc_ft{4})+numel(diff_ft{4})+numel(stat_ft{4})])]';
mot_type=[repelem([proc,diff,stat,proc,diff,stat,proc,diff,stat,proc,diff,stat], [numel(proc_ft{1}) numel(diff_ft{1}) numel(stat_ft{1}) numel(proc_ft{2}) numel(diff_ft{2}) numel(stat_ft{2}) numel(proc_ft{3}) numel(diff_ft{3}) numel(stat_ft{3}) numel(proc_ft{4}) numel(diff_ft{4}) numel(stat_ft{4})])]';

mot_typ_data_tbl=table(proc_diff_stat_all,cell_type,mot_type);

[mot_typ_p,mot_typ_tbl,mot_typ_stats]=anovan(proc_diff_stat_all,{cell_type mot_type},"varnames",["cell type","mot type"]);

[results,a,b,gnames]=multcompare(mot_typ_stats,"Dimension",[1 2],"Display","on");

mot_typ_stat_tbl = array2table(results,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
mot_typ_stat_tbl.("Group A")=gnames(mot_typ_stat_tbl.("Group A"));
mot_typ_stat_tbl.("Group B")=gnames(mot_typ_stat_tbl.("Group B"));

proc_all=[proc_ft{1}; proc_ft{2}; proc_ft{3}; proc_ft{4};];
proc_cell_type=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(proc_ft{1}) numel(proc_ft{2}) numel(proc_ft{3}) numel(proc_ft{4})])]';
proc_data_tbl=table(proc_all,proc_cell_type);

[proc_p,proc_tbl_anova,proc_stats]=anova1(proc_all,proc_cell_type,"on");

[proc_results,a,b,gnames]=multcompare(proc_stats);

proc_tbl = array2table(proc_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
proc_tbl.("Group A")=gnames(proc_tbl.("Group A"));
proc_tbl.("Group B")=gnames(proc_tbl.("Group B"));

diff_all=[diff_ft{1}; diff_ft{2}; diff_ft{3}; diff_ft{4};];
diff_cell_type=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(diff_ft{1}) numel(diff_ft{2}) numel(diff_ft{3}) numel(diff_ft{4})])]';
diff_data_tbl=table(diff_all,diff_cell_type);

[diff_p,diff_tbl_anova,diff_stats]=anova1(diff_all,diff_cell_type,"on");

[diff_results,a,b,gnames]=multcompare(diff_stats);

diff_tbl = array2table(diff_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
diff_tbl.("Group A")=gnames(diff_tbl.("Group A"));
diff_tbl.("Group B")=gnames(diff_tbl.("Group B"));

stat_all=[stat_ft{1}; stat_ft{2}; stat_ft{3}; stat_ft{4};];
stat_cell_type=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(stat_ft{1}) numel(stat_ft{2}) numel(stat_ft{3}) numel(stat_ft{4})])]';
stat_data_tbl=table(stat_all,stat_cell_type);

[stat_p,stat_tbl_anova,stat_stats]=anova1(stat_all,stat_cell_type,"on");

[stat_results,a,b,gnames]=multcompare(stat_stats);

stat_tbl = array2table(stat_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
stat_tbl.("Group A")=gnames(stat_tbl.("Group A"));
stat_tbl.("Group B")=gnames(stat_tbl.("Group B"));



%Bootstrapping for Rg
% bootstrapped_all=[];
% 
% nsamp=min([numel(rg_all{1}),numel(rg_all{2}), numel(rg_all{3}), numel(rg_all{4})]);
% 
% for dataset=1:4
% bstrp_lyso=Loic_bootstrap_code_04092019_em(rg_all{dataset},nsamp,1000,5);
% bootstrapped_dataset{dataset}=bstrp_lyso';
% end
% 
% 
% difference_30Q_45Q = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{2}.bstrap_means;
% difference_fake = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{1}.bstrap_means; %testing the bootstrapping algorithm
% difference_30Q_65Q = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{3}.bstrap_means;
% difference_30Q_81Q = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
% difference_65Q_81Q = bootstrapped_dataset{3}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
% difference_45Q_81Q = bootstrapped_dataset{2}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
% 
% figure('Name','Rg_bootstrapping_histogram','NumberTitle','off'), hold on,
% histogram(difference_30Q_45Q,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5), hold on,
% histogram(difference_30Q_65Q, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_30Q_81Q, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_65Q_81Q, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_45Q_81Q, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% % figure('Name','Rg_bootstrapping_histogram_fake','NumberTitle','off'), hold on,
% % histogram(difference_fake,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% % xlabel('Radius of Gyration (\mum)');
% % ylabel('Number of Trajectories');
% % publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_30Q_45Q','NumberTitle','off'), hold on,
% histogram(difference_30Q_45Q,'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_30Q_65Q','NumberTitle','off'), hold on,
% histogram(difference_30Q_65Q, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_30Q_81Q','NumberTitle','off'), hold on,
% histogram(difference_30Q_81Q, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_65Q_81Q','NumberTitle','off'), hold on,
% histogram(difference_65Q_81Q, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_30Q_81Q','NumberTitle','off'), hold on,
% histogram(difference_30Q_81Q, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% %Bootstrapping for alpha
% nsamp=min([numel(alph_f{1}),numel(alph_f{2}), numel(alph_f{3}), numel(alph_f{4})]);
% 
% for dataset=1:4
% bstrp_lyso_alpha=Loic_bootstrap_code_04092019_em(alph_f{dataset},nsamp,1000,5);
% bootstrapped_dataset_alpha{dataset}=bstrp_lyso_alpha';
% end
% 
% 
% difference_30Q_45Q_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{2}.bstrap_means;
% difference_fake_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{1}.bstrap_means; %testing the bootstrapping
% difference_30Q_65Q_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{3}.bstrap_means;
% difference_30Q_81Q_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{4}.bstrap_means;
% difference_65Q_81Q_alpha = bootstrapped_dataset_alpha{3}.bstrap_means - bootstrapped_dataset_alpha{4}.bstrap_means;
% difference_45Q_81Q_alpha = bootstrapped_dataset_alpha{2}.bstrap_means - bootstrapped_dataset_alpha{4}.bstrap_means;
% 
% figure('Name','Bootstrapping_histogram_alpha','NumberTitle','off'), hold on,
% histogram(difference_30Q_45Q_alpha,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5), hold on,
% histogram(difference_30Q_65Q_alpha, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_30Q_81Q_alpha, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_65Q_81Q_alpha, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_45Q_81Q_alpha, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% % figure('Name','Bootstrapping_fake_alpha','NumberTitle','off'), hold on,
% % histogram(difference_fake_alpha,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% % xlabel('\alpha');
% % ylabel('Number of Trajectories');
% % xlim([-1 1]);
% % ylim([0 2000]);
% % publication_fig(0,0,1);
% 
% 
% figure('Name','Bootstrapping_histogram_30Q_45Q_alpha','NumberTitle','off'), hold on,
% histogram(difference_30Q_45Q_alpha,'facecolor',colour_30Q, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_30Q_65Q_alpha','NumberTitle','off'), hold on,
% histogram(difference_30Q_65Q_alpha, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_30Q_81Q_alpha','NumberTitle','off'), hold on,
% histogram(difference_30Q_81Q_alpha, 'facecolor',colour_81Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_65Q_81Q_alpha','NumberTitle','off'), hold on,
% histogram(difference_65Q_81Q_alpha, 'facecolor',colour_45Q, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_45Q_81Q_alpha','NumberTitle','off'), hold on,
% histogram(difference_45Q_81Q_alpha, 'facecolor',colour_65Q, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);

%Saving all the stats and figures
tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
FolderName_1 = tempdir_1;   % Your destination folder
% %Saving stats to a table
writetable(alph_tbl,fullfile(FolderName_1, [fileprefix,'alpha.csv']),'WriteRowNames',true);
writetable(alph_data,fullfile(FolderName_1, [fileprefix,'alpha_data.csv']),'WriteRowNames',true);
writetable(rg_tbl,fullfile(FolderName_1, [fileprefix,'rg.csv']),'WriteRowNames',true);
writetable(rg_data_tbl,fullfile(FolderName_1, [fileprefix,'rg_data.csv']),'WriteRowNames',true);
writetable(dirbias_tbl,fullfile(FolderName_1, [fileprefix,'dirbias.csv']),'WriteRowNames',true);
writetable(dirbias_data,fullfile(FolderName_1, [fileprefix,'dirbias_data.csv']),'WriteRowNames',true);
writetable(mot_typ_stat_tbl,fullfile(FolderName_1, [fileprefix,'mot_typ.csv']),'WriteRowNames',true);
writetable(mot_typ_data_tbl,fullfile(FolderName_1, [fileprefix,'mot_typ_data.csv']),'WriteRowNames',true);
writetable(proc_tbl,fullfile(FolderName_1, [fileprefix,'proc.csv']),'WriteRowNames',true);
writetable(proc_data_tbl,fullfile(FolderName_1, [fileprefix,'proc_data.csv']),'WriteRowNames',true);
writetable(diff_tbl,fullfile(FolderName_1, [fileprefix,'diff.csv']),'WriteRowNames',true);
writetable(diff_data_tbl,fullfile(FolderName_1, [fileprefix,'diff_data.csv']),'WriteRowNames',true);
writetable(stat_tbl,fullfile(FolderName_1, [fileprefix,'stat.csv']),'WriteRowNames',true);
writetable(stat_data_tbl,fullfile(FolderName_1, [fileprefix,'stat_data.csv']),'WriteRowNames',true);

% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/isoHD_Neuron_motility/';  % Your destination folder
% FolderName = tempdir;  % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%  FigHandle = FigList(iFig);
%  FigName  = get(FigHandle, 'Name');
%  savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end
