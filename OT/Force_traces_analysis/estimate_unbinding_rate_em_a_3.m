% adam g. hendricks    6 oct 2021    mcgill
% Use the method described in Berger, Klumpp, and Lipowsky (2019) to
% estimate the unbinding rate (F) from stationary trap data
% 2022.02.20 - added threshold for minimum stall time
%20220403 corrected factor of 10 off in time values from previous code

clear all, %close all
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
datdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/Method_2_OT_motor_analysis/';

% force_bins = -25:0.25:25;
force_bins = -30:1:30;
Tstall_min = 0.1; %s
Fstall_min = 0; %pN
R_foff_fon_max = 0.75; %maximum ratio of off-axis to on-axis force
Fmax_fit = 20; %Maximum force to be fit for the unbinding rate exponential fit


for k_choose=1%:2
    if k_choose == 1
flist = dir(fullfile(datdir,'*WT*late*.mat'));
    elseif k_choose == 2
flist = dir(fullfile(datdir,'*S421D*late*.mat'));
    end
% 
% color_WT_plus=[0.5 0.5 0.5]; %Early
% color_WT_minus=[0 0 0]; %Early
% color_S421D_plus=[0 0.75 0.75]; %Early
% color_S421D_minus=[0 0.4 0.4]; %Early


color_WT_plus=[0.2902 0.5569  0.8863]; %Late
color_WT_minus=[0.1020    0.3490    0.6549]; %Late
color_S421D_plus=[0.4314    0.6863    0.4745]; %Late
color_S421D_minus=[0.2118    0.3843    0.2392]; %Late


kN=0; %index of Nkf
Nkf = [];
Fdetach = [];
for k = 1:numel(flist);
    datk = load(fullfile(datdir,flist(k).name));    
    % identify stall events
    for kf = 1:numel(datk.stall_ti/10),
         % calculate off-axis forces
        j_ks = find(datk.ftm/10>=datk.stall_ti(kf)/10 & datk.ftm/10<= datk.stall_tf(kf)/10); 
        [max_fpos_off_ks,j_max_fpos_off_ks] = max(abs(datk.k_trap1.*datk.fpos_off(j_ks))); %max off axis
        datk.stall_force_off(kf) = sign(datk.fpos_off(j_ks(j_max_fpos_off_ks)))2*max_fpos_off_ks;
        
        ratio_off_axis_to_on_axis{k_choose}(kf)= abs(datk.stall_force_off(kf)/datk.stall_force(kf));
        time_interval{k_choose}(kf)=datk.ftm(kf)/10;
        
        if (datk.stall_tf(kf)/10-datk.stall_ti(kf)/10) >= Tstall_min ...
                & abs(datk.stall_force(kf)) >= Fstall_min ...
                & ratio_off_axis_to_on_axis{k_choose}(kf) <= R_foff_fon_max
%                 & abs(datk.stall_force_off(kf)/datk.stall_force(kf)) <= R_foff_fon_max
            jkf = find(datk.ftm/10>datk.stall_ti(kf)/10 & datk.ftm/10<datk.stall_tf(kf)/10);
            kN = kN + 1;
            Nkf(kN,:) = hist(datk.ffm(jkf),force_bins); %count the number of time points in a given force range
            Fdetach(kN)= datk.stall_force(kf); % count the detachment events in a given force range
        end
    end 
    delta_t{k_choose}(k) = mean(diff(datk.ftm/10)); %time [s] for each data point
end
    Nkf_2{k_choose}=Nkf;
    Fdetach_2{k_choose}=Fdetach;
    Nkf=[];
    Fdetach=[];
end

N_detach_WT = hist(Fdetach_2{1},force_bins);
N_f_WT = sum(Nkf_2{1},1);
eps_WT = (1/mean(delta_t{1})).*N_detach_WT./N_f_WT;

N_detach_S421D = hist(Fdetach_2{2},force_bins);
N_f_S421D = sum(Nkf_2{2},1);
eps_S421D = (1/mean(delta_t{2})).*N_detach_S421D./N_f_S421D;

unbinding_plus_force_WT=Fdetach_2{1}(find(Fdetach_2{1}>0));
unbinding_minus_force_WT=abs(Fdetach_2{1}(find(Fdetach_2{1}<0)));
unbinding_plus_force_S421D=Fdetach_2{2}(find(Fdetach_2{2}>0));
unbinding_minus_force_S421D=abs(Fdetach_2{2}(find(Fdetach_2{2}<0)));

j_plus = find(force_bins>0);
j_minus = find(force_bins<0);

% unbinding_plus_WT=eps_WT(find(force_bins>0&N_f_WT>10&eps_WT>0&force_bins<=Fmax_fit)); % Setting a threshold for the number of detachment events
% unbinding_minus_WT=eps_WT(find(force_bins<0&N_f_WT>10&eps_WT>0&abs(force_bins)<=Fmax_fit));
% unbinding_plus_S421D=eps_S421D(find(force_bins>0&N_f_S421D>10&eps_S421D>0&force_bins<=Fmax_fit)); 
% unbinding_minus_S421D=eps_S421D(find(force_bins<0&N_f_S421D>10&eps_S421D>0&abs(force_bins)<=Fmax_fit));
% plus_force_bins_WT=force_bins(find(force_bins>0&N_f_WT>10&eps_WT>0&force_bins<=Fmax_fit)).';
% minus_force_bins_WT=force_bins(find(force_bins<0&N_f_WT>10&eps_WT>0&abs(force_bins)<=Fmax_fit)).';
% plus_force_bins_S421D=force_bins(find(force_bins>0&N_f_S421D>10&eps_S421D>0&force_bins<=Fmax_fit)).';
% minus_force_bins_S421D=force_bins(find(force_bins<0&N_f_S421D>10&eps_S421D>0&abs(force_bins)<=Fmax_fit)).';

% 
unbinding_plus_WT=eps_WT(find(force_bins>0&N_f_WT>10&eps_WT>0&eps_WT<50&force_bins<Fmax_fit)); % Setting a threshold for the number of detachment events
unbinding_minus_WT=eps_WT(find(force_bins<0&N_f_WT>10&eps_WT>0&eps_WT<50&abs(force_bins)<Fmax_fit));
unbinding_plus_S421D=eps_S421D(find(force_bins>0&N_f_S421D>10&eps_S421D>0&eps_S421D<50&force_bins<Fmax_fit)); 
unbinding_minus_S421D=eps_S421D(find(force_bins<0&N_f_S421D>10&eps_S421D>0&eps_S421D<50&abs(force_bins)<Fmax_fit));
plus_force_bins_WT=force_bins(find(force_bins>0&N_f_WT>10&eps_WT>0&eps_WT<50&force_bins<Fmax_fit)).';
minus_force_bins_WT=force_bins(find(force_bins<0&N_f_WT>10&eps_WT>0&eps_WT<50&abs(force_bins)<Fmax_fit)).';
plus_force_bins_S421D=force_bins(find(force_bins>0&N_f_S421D>10&eps_S421D>0&eps_S421D<50&force_bins<Fmax_fit)).';
minus_force_bins_S421D=force_bins(find(force_bins<0&N_f_S421D>10&eps_S421D>0&eps_S421D<50&abs(force_bins)<Fmax_fit)).';

% fit_WT_plus=fit(plus_force_bins_WT(1:end-1),unbinding_plus_WT(1:end-1).','exp1');
% fit_S421D_plus=fit(plus_force_bins_S421D(1:end-1),unbinding_plus_S421D(1:end-1).','exp1');
% fit_WT_minus=fit(minus_force_bins_WT(2:end),unbinding_minus_WT(2:end).','exp1');
% fit_S421D_minus=fit(minus_force_bins_S421D(2:end),unbinding_minus_S421D(2:end).','exp1');

fit_WT_plus=fit(plus_force_bins_WT,unbinding_plus_WT.','exp1');
fit_S421D_plus=fit(plus_force_bins_S421D,unbinding_plus_S421D.','exp1');
fit_WT_minus=fit(minus_force_bins_WT,unbinding_minus_WT.','exp1');
fit_S421D_minus=fit(minus_force_bins_S421D,unbinding_minus_S421D.','exp1');

figure('Name','Unbinding_rate_force_dependence','NumberTitle','off'), hold on, 
plot(force_bins,eps_WT,'o','linewidth',2,'Color',color_WT_plus),hold on
plot(force_bins,eps_S421D,'o','linewidth',2,'Color',color_S421D_plus);
xlabel('Force (pN)'), ylabel('Unbinding Rate (1/s)')
set(gca,'yscale','log');
publication_fig(0,0,1);

figure('Name','Unbinding_rate_force_dependence 2','NumberTitle','off'), hold on, 
plot(plus_force_bins_WT,unbinding_plus_WT,'o','linewidth',2,'Color',color_WT_plus),hold on
plot(minus_force_bins_WT,unbinding_minus_WT,'o','linewidth',2,'Color',color_WT_plus),hold on
plot(plus_force_bins_S421D,unbinding_plus_S421D,'o','linewidth',2,'Color',color_S421D_plus), hold on
plot(minus_force_bins_S421D,unbinding_minus_S421D,'o','linewidth',2,'Color',color_S421D_plus), hold on
WT1=plot(fit_WT_plus,plus_force_bins_WT,unbinding_plus_WT);
S421D1=plot(fit_S421D_plus,plus_force_bins_S421D,unbinding_plus_S421D);
WT2=plot(fit_WT_minus,minus_force_bins_WT,unbinding_minus_WT);
S421D2=plot(fit_S421D_minus,minus_force_bins_S421D, unbinding_minus_S421D);
set([WT1 WT2], 'Color',color_WT_plus);
set([S421D1 S421D2], 'Color',color_S421D_plus);
xlabel('Force (pN)'), ylabel('Unbinding Rate (1/s)')
xlim([-20 20]);
% set(gca,'yscale','log');
legend('off');
publication_fig(0,0,1);

figure('Name','Unbinding_rate_force_dependence 3','NumberTitle','off'), hold on, 
plot(plus_force_bins_WT,unbinding_plus_WT,'o','linewidth',2,'Color',color_WT_plus),hold on
plot(minus_force_bins_WT,unbinding_minus_WT,'o','linewidth',2,'Color',color_WT_plus),hold on
plot(plus_force_bins_S421D,unbinding_plus_S421D,'o','linewidth',2,'Color',color_S421D_plus), hold on
plot(minus_force_bins_S421D,unbinding_minus_S421D,'o','linewidth',2,'Color',color_S421D_plus), hold on
WT1=plot(fit_WT_plus);
S421D1=plot(fit_S421D_plus);
WT2=plot(fit_WT_minus);
S421D2=plot(fit_S421D_minus);
set([WT1 WT2], 'Color',color_WT_plus);
set([S421D1 S421D2], 'Color',color_S421D_plus);
xlabel('Force (pN)'), ylabel('Unbinding Rate (1/s)')
xlim([-20 20]);
% set(gca,'yscale','log');
legend('off');
publication_fig(0,0,1);

figure('Name','Force_distribution_unbinding','NumberTitle','off'), hold on, 
b1=bar(force_bins(j_plus),N_detach_WT(j_plus)./sum(N_detach_WT),1);
b1.FaceAlpha=0.3;
b1.FaceColor=color_WT_minus;
b2=bar(-force_bins(j_minus), -N_detach_WT(j_minus)./sum(N_detach_WT),1);
b2.FaceAlpha=0.3;
b2.FaceColor=color_WT_minus;
b3=bar(force_bins(j_plus),N_detach_S421D(j_plus)./sum(N_detach_S421D),1);
b3.FaceAlpha=0.3;
b3.FaceColor=color_S421D_minus;
b4=bar(-force_bins(j_minus), -N_detach_S421D(j_minus)./sum(N_detach_S421D),1);
b4.FaceAlpha=0.3;
b4.FaceColor=color_S421D_minus;
xlabel('Force (pN)')
ylabel('Fraction')
publication_fig(0,0,1);

%Emily adding the error bars to the force distribution plot

mean_ant_WT=mean(unbinding_plus_force_WT);
mean_ant_S421D=mean(unbinding_plus_force_S421D);
mean_ret_WT=mean(unbinding_minus_force_WT);
mean_ret_S421D=mean(unbinding_minus_force_S421D);

sd_ant_WT=std(unbinding_plus_force_WT);
sd_ant_S421D=std(unbinding_plus_force_S421D);
sd_ret_WT=std(unbinding_minus_force_WT);
sd_ret_S421D=std(unbinding_minus_force_S421D);

sem_ant_WT=std(unbinding_plus_force_WT)/sqrt(numel(unbinding_plus_force_WT));
sem_ant_S421D=std(unbinding_plus_force_S421D)/sqrt(numel(unbinding_plus_force_S421D));
sem_ret_WT=std(unbinding_minus_force_WT)/sqrt(numel(unbinding_minus_force_WT));
sem_ret_S421D=std(unbinding_minus_force_S421D)/sqrt(numel(unbinding_minus_force_S421D));

%Force histograms for all stalls above the threshold 0.25s stall force,
%binned by 0.5pN bins 
figure('Name','Force histogram both conditions','NumberTitle','off'), hold on,
xlabel('Force (pN)'), ylabel('Fraction of Events'),
hold on, 
b1=bar(force_bins(j_plus),N_detach_S421D(j_plus)./sum(N_detach_S421D),1,'facecolor',color_S421D_minus,'edgecolor','k','linewidth',0.5);
b1.FaceAlpha=0.3;
hold on, 
b2=bar(-force_bins(j_minus), -N_detach_S421D(j_minus)./sum(N_detach_S421D),1,'facecolor',color_S421D_minus,'edgecolor','k','linewidth',0.5);
b2.FaceAlpha=0.3;
hold on,
b3=bar(force_bins(j_plus),N_detach_WT(j_plus)./sum(N_detach_WT),1,'facecolor',color_WT_plus,'edgecolor','k','linewidth',0.5);
b3.FaceAlpha=0.3;
hold on, 
b4=bar(-force_bins(j_minus), -N_detach_WT(j_minus)./sum(N_detach_WT),1,'facecolor',color_WT_plus,'edgecolor','k','linewidth',0.5);
b4.FaceAlpha=0.3;
hold on, scatter(mean_ant_WT,0.15,100,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor','none');
hold on, plot([mean_ant_WT-sd_ant_WT; mean_ant_WT+sd_ant_WT], [0.15; 0.15],'Color', color_WT_minus);
hold on, scatter(mean_ret_WT,-0.15,100,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor','none');
hold on, plot([mean_ret_WT-sd_ret_WT; mean_ret_WT+sd_ret_WT], [-0.15; -0.15],'Color', color_WT_minus);
hold on, scatter(mean_ant_S421D,0.2,100,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor','none');
hold on, plot([mean_ant_S421D-sd_ant_S421D; mean_ant_S421D+sd_ant_S421D], [0.2; 0.2], 'Color', color_S421D_plus);
hold on, scatter(mean_ret_S421D,-0.2,100,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor','none');
hold on, plot([mean_ret_S421D-sd_ret_S421D; mean_ret_S421D+sd_ret_S421D], [-0.2; -0.2],'Color', color_S421D_plus);
xlim([-1.5 25]);
ylim([-0.25 0.25]);
publication_fig(0,0,1);
% axis tight

population_plus_WT=N_detach_WT(j_plus)./sum(N_detach_WT);
population_minus_WT=N_detach_WT(j_minus)./sum(N_detach_WT);
population_plus_S421D=N_detach_S421D(j_plus)./sum(N_detach_S421D);
population_minus_S421D=N_detach_S421D(j_minus)./sum(N_detach_S421D);
% [kstest_plus_minus,kstest_plus_minus_p]=kstest2(population_plus,population_minus)

plus_forces_fraction_WT=sum(N_detach_WT(j_plus)./sum(N_detach_WT));
minus_forces_fraction_WT=sum(-N_detach_WT(j_minus)./sum(N_detach_WT));
plus_forces_fraction_S421D=sum(N_detach_S421D(j_plus)./sum(N_detach_S421D));
minus_forces_fraction_S421D=sum(-N_detach_S421D(j_minus)./sum(N_detach_S421D));

[WT_plus_vs_minus_f, p_WT_plus_vs_minus_f]=kstest2(unbinding_plus_force_WT,abs(unbinding_minus_force_WT))
[S421D_plus_vs_minus_f, p_S421D_plus_vs_minus_f]=kstest2(unbinding_plus_force_S421D,abs(unbinding_minus_force_S421D))
[WT_S421D_plus_u_f, p_WT_S421D_plus_u_f]=kstest2(unbinding_plus_force_WT,unbinding_plus_force_S421D)
[WT_S421D_minus_u_f, p_WT_S421D_minus_u_f]=kstest2(abs(unbinding_minus_force_WT),abs(unbinding_minus_force_S421D))

WT_plus_bstrp=Loic_bootstrap_code_04092019_em(unbinding_plus_force_WT, numel(unbinding_plus_force_WT),1000,0.05)
S421D_plus_bstrp=Loic_bootstrap_code_04092019_em(unbinding_plus_force_S421D, numel(unbinding_plus_force_S421D),1000,0.05)
WT_minus_bstrp=Loic_bootstrap_code_04092019_em(unbinding_minus_force_WT, numel(unbinding_minus_force_WT),1000,0.05)
S421D_minus_bstrp=Loic_bootstrap_code_04092019_em(unbinding_minus_force_S421D, numel(unbinding_minus_force_S421D),1000,0.05)

%plotting bootstrap difference histograms
figure, histogram(S421D_plus_bstrp.bstrap_means-WT_plus_bstrp.bstrap_means);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('S421D-WT Plus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(WT_plus_bstrp.bstrap_means-abs(WT_minus_bstrp.bstrap_means));
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('WT Plus-WT Minus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(abs(S421D_minus_bstrp.bstrap_means)-abs(WT_minus_bstrp.bstrap_means));
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('S421D-WT Minus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(S421D_plus_bstrp.bstrap_means-abs(S421D_minus_bstrp.bstrap_means));
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('S421D Plus-S421D Minus'), ylabel('Frequency');
publication_fig(0,0,1);

ci_bstrp_WT_plus=quantile(WT_plus_bstrp.bstrap_means, [0.05 0.95]);
ci_bstrp_S421D_plus=quantile(S421D_plus_bstrp.bstrap_means, [0.05 0.95]);
ci_bstrp_WT_minus=quantile(WT_minus_bstrp.bstrap_means, [0.05 0.95]);
ci_bstrp_S421D_minus=quantile(S421D_minus_bstrp.bstrap_means, [0.05 0.95]);

p_WT_S421D_plus=numel(find(WT_plus_bstrp.bstrap_means-S421D_plus_bstrp.bstrap_means>0))/numel(WT_plus_bstrp.bstrap_means-S421D_plus_bstrp.bstrap_means)
p_WT_plus_minus=numel(find(WT_plus_bstrp.bstrap_means-abs(WT_minus_bstrp.bstrap_means)>0))/numel(WT_plus_bstrp.bstrap_means-abs(WT_minus_bstrp.bstrap_means))
p_WT_S421D_minus=numel(find(abs(WT_minus_bstrp.bstrap_means)-abs(S421D_minus_bstrp.bstrap_means)>0))/numel(abs(WT_minus_bstrp.bstrap_means)-abs(S421D_minus_bstrp.bstrap_means))
p_S421D_plus_minus=numel(find(S421D_plus_bstrp.bstrap_means-abs(S421D_minus_bstrp.bstrap_means)>0))/numel(S421D_plus_bstrp.bstrap_means-abs(S421D_minus_bstrp.bstrap_means))

figure('Name','Force histogram both conditions 2','NumberTitle','off'), hold on,
xlabel('Force (pN)'), ylabel('Fraction of Events'),
hold on, 
b1=bar(force_bins(j_plus),N_detach_S421D(j_plus)./sum(N_detach_S421D),1,'facecolor',color_S421D_minus,'edgecolor','k','linewidth',0.5);
b1.FaceAlpha=0.3;
hold on, 
b2=bar(-force_bins(j_minus), -N_detach_S421D(j_minus)./sum(N_detach_S421D),1,'facecolor',color_S421D_minus,'edgecolor','k','linewidth',0.5);
b2.FaceAlpha=0.3;
hold on,
b3=bar(force_bins(j_plus),N_detach_WT(j_plus)./sum(N_detach_WT),1,'facecolor',color_WT_plus,'edgecolor','k','linewidth',0.5);
b3.FaceAlpha=0.3;
hold on, 
b4=bar(-force_bins(j_minus), -N_detach_WT(j_minus)./sum(N_detach_WT),1,'facecolor',color_WT_plus,'edgecolor','k','linewidth',0.5);
b4.FaceAlpha=0.3;
hold on, scatter(WT_plus_bstrp.mean_of_means,0.15,50,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor','none');
hold on, plot([ci_bstrp_WT_plus(1); ci_bstrp_WT_plus(2)], [0.15; 0.15],'Color', color_WT_minus);
hold on, scatter(WT_minus_bstrp.mean_of_means,-0.15,50,'MarkerFaceColor',color_WT_minus,'MarkerEdgeColor','none');
hold on, plot([ci_bstrp_WT_minus(1); ci_bstrp_WT_minus(2)], [-0.15; -0.15],'Color', color_WT_minus);
hold on, scatter(S421D_plus_bstrp.mean_of_means,0.2,50,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor','none');
hold on, plot([ci_bstrp_S421D_plus(1); ci_bstrp_S421D_plus(2)], [0.2; 0.2], 'Color', color_S421D_plus);
hold on, scatter(S421D_minus_bstrp.mean_of_means,-0.2,50,'MarkerFaceColor',color_S421D_plus,'MarkerEdgeColor','none');
hold on, plot([ci_bstrp_S421D_minus(1); ci_bstrp_S421D_minus(2)], [-0.2; -0.2],'Color', color_S421D_plus);
xlim([-1.5 25]);
ylim([-0.25 0.25]);
publication_fig(0,0,1);

% fileprefix='20220402_late_stall_0pt1_unbinding_';
% % 
% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end
