%% Plot x,y per cell and calculate tracking uncertainty
%% Emily N. P. Prowse - 2024/01/20 based on Abdullah R. Chaudhary - 05/16/2020

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');

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
htraj=gcf;
cmap=colormap(parula(1000));
dmin=0;
stat_disp=[];
%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for k_choose = 1:6
    
if k_choose == 1    % Condition 1
    col1=[0.5 0.5 0.5];
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    nm='WT';
    
elseif k_choose == 2   % Condition 2
    col1=[0 0.4 0.4];
    cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_lyso_mats');
    fl='*.mat';
fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 3
    cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_lyso_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 4
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 5
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_lyso_IFg_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 6
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_lyso_IFg_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 7
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_mats/');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 8
    cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_bdnf_mats/');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 9
    cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_bdnf_mats/');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose ==10
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_mats/');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 11
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_bdnf_IFg_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 12
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_bdnf_IFg_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 13
    cd('/Volumes/Emily_htt_2/Neuron/isoHD30Q/30Q_mito_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 14
    cd('/Volumes/Emily_htt_2/Neuron/isoHD45Q/45Q_mito_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 15
    cd('/Volumes/Emily_htt_2/Neuron/isoHD65Q/65Q_mito_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
elseif k_choose == 16
    cd('/Volumes/Emily_htt_2/Neuron/isoHD81Q/81Q_mito_mats');
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;    
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
    for j=1:numel(ab) %j is the trajectory number
%         numel(ab(j).position)
        if numel(ab(j).position)>15 %Emily changed this value from 75 to 0 on 15/08/2020
        kf=kf+1;
        xk{kf}=ab(j).xk;
        yk{kf}=ab(j).yk;
        tk{kf}=ab(j).tk;
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
        %Calculating track displacement
         displacement{kf} = sqrt((xk{kf}-xk{kf}(1)).^2 + (yk{kf}-yk{kf}(1)).^2);
         track_displacement(kf) = displacement{kf}(end);%distance between beginning and end point (um)
         pdk = polyfit(tk{kf},displacement{kf},1);
         track_velocity(kf) = pdk(1);
        
        if (0<alpha_1)<0.5 && track_displacement(kf)<0.5  
            stat_disps{kf}=track_displacement(kf);
            rg_0=func_Rg_Linda_v2(xk{kf},yk{kf});
            rg_1=[rg_1,rg_0];
            dat_pos_tim{kf}=smooth(ab(j).position,30,'sgolay',1);
            xk_yk{kf}=[xk{kf},yk{kf}];
           
%             figure(k), hold on, plot(xk{kf},yk{kf},'.-','Color',cmap(kf,:)); %plot all
        else
        rg_1=[];
        dat_pos_tim{kf}=[];
        xk_yk{kf}=[];
        stat_disps{kf}=[];
        end
        
        else
        end
      
    end

    stat_disp{k}=cell2mat(stat_disps(~cellfun('isempty',stat_disps)));
    all_disp{k}=track_displacement;
%     pos=[pos,dat_pos_tim];%commmented out Jan 20 2024
%     position{k}=pos;%commmented out Jan 20 2024
    
    %figure(k), hold on, 
    %cellfun(@(x,y) plot(x,y, 'Color','k'), xk_yk);
%     fig=figure(k); %Em commented out to try another style of figure

%     fig=figure('Name',fln,'NumberTitle','off');
%     ax = gca(fig);
% %     ax.ColorOrder=[0.7 0.7 0.7; 0.7 0.7 0.7;]; %WT r5
%  ax.ColorOrder=[0.6 0.87 0.87; 0.6 0.87 0.87];%S421D r5
% %     ax.ColorOrder=[0.729 0.831 0.957; 0.35 0.4 0.45; 0.549 0.667 0.745]; %WT lyso
% %  ax.ColorOrder=[0.729 0.831 0.957; 0.729 0.831 0.957]; %WT lyso
% %     ax.ColorOrder=[0.757 0.867 0.776;0.757 0.867 0.776];%S421D lyso
% %     ax.ColorOrder=[0 0.9 0.9; 0 0.4 0.4; 0 0.75 0.75; 0 0.55 0.55];%S421D r5
% %     ax.ColorOrder=[0.729 0.831 0.957; 0.35 0.4 0.45; 0.549 0.667 0.745]; %WT lyso
% %     ax.ColorOrder=[0.757 0.867 0.776; 0.35 0.4 0.35; 0.522 0.74 0.549];%S421D lyso
% 
%     hold on
%     cellfun(@(X,Y) plot(ax, X, Y), xk,yk) %Uncomment to plot xy positions of trajectories
% %     cellfun(@(X,Y) plot(ax, X, Y), tim, pos_tim) %Uncomment to plot
% %     position over time
%     publication_fig(0,0,1);
% %     xlabel('xk (\mum)'); ylabel('yk (\mum)');
%     xlabel('Time (s)'); ylabel('Position (\mum)');
% %     axis equal
%     hold off
    

%     fl_nm=[fls(k).name, nm];
%     saveas(gca,fl_nm,'png')
    
    clear res  xk yk tk rg_2 xk_yk tim pos_tim; %Emily added tim and pos_tim Mar 28/22, dat_pos_tim
    kf=0;

    


end

stat=cell2mat(stat_disp);
all=cell2mat(all_disp);

figure(k_choose)
xdisp = 0:0.1:5;
ndisp = hist(all,xdisp);
jk = find(all>=dmin);
[phat, pci] = mle(all(jk)-dmin,'distribution','exp');
xx = 0:0.1:xdisp(end);
yy = exppdf(xx,phat);

figure(k_choose), bar(xdisp, ndisp./sum(ndisp))
hold on, plot(xx,1.7.*yy./sum(yy),'r-','linewidth',2)
xlabel('Track Displacement (\mum)'), ylabel('Fraction of Runs')
%gtext(['Mean = ',num2str(phat,3),' \mum'])

stat_per_cond{k_choose}=stat;
all_per_cond{k_choose}=all;

clear stat_disp all_disp
    
end

stat_all=cell2mat(stat_per_cond);
disp_all=cell2mat(all_per_cond);

figure(100),
xdisp = 0:0.1:5;
ndisp = hist(disp_all,xdisp);
jk = find(disp_all>=dmin);
[phat, pci] = mle(disp_all(jk)-dmin,'distribution','exp');
xx = 0:0.1:xdisp(end);
yy = exppdf(xx,phat);

figure(100), bar(xdisp, ndisp./sum(ndisp))
hold on, plot(xx,1.7.*yy./sum(yy),'r-','linewidth',2)
xlabel('Track Displacement (\mum)'), ylabel('Fraction of Runs')
pbaspect([2 1 1]);
%gtext(['Mean = ',num2str(phat,3),' \mum'])

options = statset('Display','final','MaxIter',1000,'TolX',1e-4,'robust','on'); %options for the gaussian mixture distribution fit, changed from 1e-8 to 1e-4


xdisp = 0:0.025:0.5;
ndisp = hist(stat_all,xdisp);
jk = find(stat_all>=dmin);
[m,s,mci,sci]=normfit(stat_all',0.05);
xx = 0:0.1:xdisp(end);
yy = normpdf(xx,m,s);

figure(200), bar(xdisp, ndisp./sum(ndisp));
% hold on, plot(xx,0.6.*yy./(sum(yy)),'r-','linewidth',2);
xlabel('Track Displacement (\mum)'), ylabel('Fraction of Runs');

for kdist = 4 %The number of gaussian fits to perform, this can be changed to a range and optimized by minimizing the BIC criteria
    fit = gmdistribution.fit(stat_all',kdist,'Options',options,'Start','randSample'); %Fitting the gaussian distributions
    bic(kdist)=fit.BIC;
    pdf_fit=pdf(fit,xdisp');  %Probability density function for the gaussian fit (for plotting)
    figure(1000)
    plot(bic,'b--');
    figure(10000)
    plot(xdisp,0.025*pdf_fit,'Color','r','LineStyle','-','LineWidth', 3);
end

clear ab xk yk pos kf;


toc
