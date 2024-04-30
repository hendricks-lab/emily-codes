%% Abdullah R. Chaudhary
%% This code will only run .csv files
clear all; close all; clc; 
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
set(0,'DefaultFigureWindowStyle','docked');

tic
% for k_choose = 1:2
for k_choose = 4

if k_choose == 2    % 30Q mito
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_tracks/');
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_mats/';
elseif k_choose == 3    % 45Q mito
     cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_tracks/');
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_mats/';
elseif k_choose == 4    % 65Q mito
     cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_tracks/');
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_mats/';
elseif k_choose == 5    % 81Q mito
     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_tracks/');
    save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_mats/';
end
% 
% if k_choose == 1    % 18Q bdnf
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_bdnf_mats/';
%     cd('/Volumes/Emily_2022/');
% %     save_dir='/Volumes/Emily_2022/Fake_TM_flux_test_mat/';
% elseif k_choose == 2    % 30Q bdnf
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_mats/';
%  elseif k_choose == 3    % 45Q bdnf
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_BDNF_mats/';
% elseif k_choose == 4    % 65Q bdnf
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_mats/';
% elseif k_choose == 5   % 81Q bdnf
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_mats/';
% end
% 
% if k_choose == 1 % 30Q bdnf +IFg
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_mats/';
% elseif k_choose == 2   % 81Q bdnf+IFg
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_mats/';
% elseif k_choose == 3   % 30Q lyso+IFg
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_mats/';
% elseif k_choose == 4   % 81Q lyso+IFg
%      cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_mats/';
% end

% if k_choose == 1    % 18Q lyso
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD18Q/18Q_lyso_mats/';
% elseif k_choose == 2    % 30Q lyso
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_mats/';
% elseif k_choose == 3   % 45Q lyso
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_mats/';
% elseif k_choose == 4    % 65Q lyso
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_mats/';
% elseif k_choose == 5    % 81Q lyso
%     cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_tracks/');
%     save_dir='/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats/';
% end

% Asterik indicates all the files with names starting with Spots and with
% extension .csv
Mot_file='Spots*.csv';
DT=0.12;    % single channel exposure time 120ms
%% Variables defined
tm_bw_rev=[];
rn_bw_rev=[];
dat=[];

dmot=dir(Mot_file); % Pick the trajectory files in .mat format

for k=1:length(dmot) % k is a vector of length dmot (no k=4)
    display(k)
dat=readmatrix(dmot(k).name); %load everything in dmot
[filepath,name,ext] = fileparts(dmot(k).name);
display(dmot(k).name);  % display
cmap = colormap(lines(100));
Ndat=max(dat(1:end,3))+1; %Emily modifying for mitochondrial tracking analysis because some tracks don't seem to have trackIDs in the third column (20221003)
kp=0; 

for kd=1:Ndat
    kp=kp+1;
    jkd=find(dat(:,3)==(kd-1)); %indices of track_k

    xk=dat(jkd,5); %x position (um)
    yk=dat(jkd,6); %y position (um)
    timek=[1:numel(xk)].*DT;    % Time (s)
%     origin_x=dat(2,22); % Point of origin in x-axis old trackmate 
%     origin_y=dat(2,23); % Point of origin in y-axis old trackmate 
%     axon_length=dat(2,24); %Measured axon length in um.
    %For mitochondrial analysis
    origin_x=dat(5,32); % Point of origin in x-axis weka segmentation trackmate -mito data
    origin_y=dat(5,33); % Point of origin in y-axis weka segmentation trackmate -mito data
    axon_length=dat(5,34); %Measured axon length in um. 
    mito_area=dat(jkd,27); %Mitochondrial area, comment out for other cargoes
    
    cargo_intensity=dat(jkd,17); %For all

%     origin_x=dat(5,21); % Point of origin in x-axis -BDNF/lyso data (new trackmate)
%     origin_y=dat(5,22); % Point of origin in y-axis -BDNF/lyso data (new trackmate)
%     axon_length=dat(5,23); % Measured axon length in um.
%     
    %if range(xk)>=0.5 || range(yk)>=0.5
    % Determine median position from kcluster analysis:
    [idx,kclust]=kmeans([xk,yk],1);
    dx_med=(origin_x-kclust(1));
    dy_med=(origin_y-kclust(2));
    
    % Vectors of x and y components:
    dx=[0;diff(xk)]; dy=[0;diff(yk)];
    dx_origin=(origin_x-xk); dy_origin=(origin_y-yk);
    
    for j=1:numel(dx)
    dt_prd(j)=dot([dx(j),dy(j)],[dx_origin(j),dy_origin(j)],2);
    end
    
    mag_vect=sqrt(dx_origin.^2 + dy_origin.^2);
    dt_prd_norm=dt_prd'./mag_vect;
    position=cumsum(dt_prd_norm);

    %Saving the variables
    ab(kp).position=position.*-1;
    ab(kp).xk=xk;
    ab(kp).yk=yk;
    ab(kp).tk=timek;
    ab(kp).med_clust_x=dx_med;
    ab(kp).med_clust_y=dy_med;
    ab(kp).origin_x=origin_x;
    ab(kp).origin_y=origin_y;
    ab(kp).axon_length=axon_length;
    ab(kp).mito_area=mito_area; %Comment out for other cargos (only for mito)
    ab(kp).cargo_intensity=cargo_intensity;
    clear dt_prd;
    %else
%     ab(kp).position=[];
%     ab(kp).xk=[];
%     ab(kp).yk=[];
%     ab(kp).tk=[];
%     ab(kp).med_clust_x=[];
%     ab(kp).med_clust_y=[];
%     ab(kp).origin_x=[];
%     ab(kp).origin_y=[];
    %end
  
end
    save([save_dir, num2str(k_choose) name],'ab');
    
    clear ab;
end


end
toc


