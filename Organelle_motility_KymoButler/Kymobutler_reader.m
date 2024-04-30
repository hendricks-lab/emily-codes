clear all;
tic
addpath('/Volumes/Emily_htt_2/New_General_codes/');

colour_30Q=[0.5 0.5 0.5];
colour_45Q=[0.75 0 0];
colour_65Q=[0.5 0 0];
colour_81Q=[0.25 0 0];

DT=0.12;    % single channel exposure time 120ms

%% Variables defined
for k_choose = 1:4

if k_choose == 1    % 30Q
    %30Q BDNF
%     cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_bdnf_mats/';
%30Q mitochondria
    cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_tracks/');
    save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_mito_mats/';
    %30Q Lysosomes
%     cd('/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_KB/30Q_KB_lyso_mats/';

elseif k_choose == 2    
    % 45Q BDNF
%      cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_bdnf_mats/';
    %45Q mitochondria
    cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_tracks/');
    save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_mito_mats/';
%45Q lysosomes
%     cd('/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/45Q_KB/45Q_KB_lyso_mats/';

elseif k_choose == 3   % 65Q
    %65Q BDNF
%      cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_bdnf_mats/';
    %65Q mitochondria
    cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_tracks/');
    save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_mito_mats/';
    %65Q lyso
%     cd('/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/65Q_KB/65Q_KB_lyso_mats/';
elseif k_choose == 4    % 81Q
    %81Q BDNF
%      cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_bdnf_mats/';
    %81Q mitochondria
    cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_tracks/');
    save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_mito_mats/';
    %81Q lysosome
%     cd('/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_KB/81Q_KB_lyso_mats/';
elseif k_choose == 5 
        %30Q BDNF IFg
%      cd('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_bdnf_mats/';
%       30Q Lyso IFg
%      cd('/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/30Q_IFg_KB/30Q_IFg_KB_lyso_mats/';
elseif k_choose== 6 
    %81Q BDNF IFg
%      cd('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_bdnf_mats/';
    %81Q lyso IFg
%      cd('/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_tracks/');
%     save_dir='/Volumes/Emily_2022/Neurons/81Q_IFg_KB/81Q_IFg_KB_lyso_mats/';
end

Mot_file='trackCoordinates*.csv';
DT=0.12;    % single channel exposure time 120ms
scale=1; %0.16; %this should be incorporated already

dat=[];

dmot=dir(Mot_file); % Pick the trajectory files in .mat format

for k=1:length(dmot) % k is a vector of length dmot
    disp(k)
    dat=readmatrix(dmot(k).name); %load everything in dmot
    [filepath,name,ext] = fileparts(dmot(k).name);
    display(dmot(k).name);  % display
    cmap = colormap(lines(100000));
    Ndat=size(dat,2); %need the number of columns from the csv
    kp=0; 
    
for kd=1:Ndat/2
    kp=kp+1;

    posk=dat(:,kd*2).*scale; % Position (μm)
    timek=(1:(numel(dat(:,kd*2-1)))).*DT; %Time (s)

    pos=posk(~isnan(posk)); % Position (μm), need to remove NaNs
    time=timek(1:numel(pos)); %Time (s), need to remove NaNs
    realtimek=dat(:,kd*2-1); %Time from raw data (not always starting at zero)
    realtime=realtimek(~isnan(realtimek));

    kbpos(kp).position=pos;
    kbpos(kp).tk=time';
    kbpos(kp).rt=realtime;

%     figure(k*10) %Plotting the kymograph
%     plot(pos,time,'-','Color',cmap(kd,:),'linewidth',2)
%     xlabel('Position (\mum)'), ylabel('Time (s)'), hold on
%     set(gca,'Ydir','reverse')
end
    save([save_dir, num2str(k_choose) name],'kbpos');
    
    clear kbpos;
end

end
toc