% adam g. hendricks    6 oct 2021    mcgill
% Use the method described in Berger, Klumpp, and Lipowsky (2019) to
% estimate the unbinding rate (F) from stationary trap data
% 2022.02.20 - added threshold for minimum stall time

clear all, %close all
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/'); %reference the folder containing the fit_single_exponential function

%datdir = '/Users/adhendri/OneDrive - McGill University/hendricks-lab-share/data/EP_force_traces/ctl_fig_trace/'
%datdir = '/Users/adhendri/OneDrive - McGill University/hendricks-lab-share/data/EP_force_traces/tau_fig_trace/'
%time_scale_factor = 1;

% datdir = '/Users/adhendri/OneDrive - McGill University/hendricks-lab-share/data/Emily OT/Full_length_force_trace_data/S421D_early/'
% datdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/Method_2_OT_motor_analysis/';
datdir = '/Volumes/Emily_2022/Omar_OT/U2OS/Good_force_traces/Force_mat/'
% 'WT_early', 'S421D_early', 'WT_late', 'S421D_late')
time_scale_factor = 0.1 %corrects time scale in Emily's data (Fs = 20 kHz)
flist = dir(fullfile(datdir,'*.mat'));
% color = [0.5 0.5 0.5] %gray - WT early
% alpha = 1;
% color = [0, 0.75, 0.75]; %green - S421D early
% alpha = 0.5;
% color = [0.2902 0.5569  0.8863] %blue - WT late
% alpha = 1;
color = [0.4314    0.6863    0.4745] %green - S421D late
alpha = 0.5;



force_bins = -30:1:30;
Tstall_min = 0.1; %s
Fstall_min = 0; %pN
R_foff_fon_max = 0.65; %maximum ratio of off-axis to on-axis force
N_f_min = 10; %minumum number of time points to include in unbinding rate calculation
Fmax_fit = 20;%max(force_bins); %maximum force to include in exponential fit

kN=0; %index of Nkf
Nkf = [];
Fdetach = [];
for k = 1:numel(flist);
    datk = load(fullfile(datdir,flist(k).name));
    datk.stall_ti = datk.stall_ti.*time_scale_factor;
    datk.stall_tf = datk.stall_tf.*time_scale_factor;
    datk.ftm = datk.ftm.*time_scale_factor;
    datk.stall_time = datk.stall_time.*time_scale_factor;
    datk.stall_vel = datk.stall_vel./time_scale_factor;
    
    % identify stall events
    for kf = 1:numel(datk.stall_ti),
         % calculate off-axis forces
        j_ks = find(datk.ftm>=datk.stall_ti(kf) & datk.ftm<= datk.stall_tf(kf)); 
        [max_fpos_off_ks,j_max_fpos_off_ks] = max(abs(datk.k_trap1.*datk.fpos_off(j_ks)));
        datk.stall_force_off(kf) = sign(datk.fpos_off(j_ks(j_max_fpos_off_ks)))*max_fpos_off_ks;
        
        if datk.stall_tf(kf)-datk.stall_ti(kf) >= Tstall_min ...
                & abs(datk.stall_force(kf)) >= Fstall_min ...
                & abs(datk.stall_force_off(kf)/datk.stall_force(kf)) <= R_foff_fon_max
            jkf = find(datk.ftm>datk.stall_ti(kf) & datk.ftm<datk.stall_tf(kf));
            kN = kN + 1;
            Nkf(kN,:) = hist(datk.ffm(jkf),force_bins); %count the number of time points in a given force range
            Fdetach(kN) = datk.stall_force(kf); % count the detachment events in a given force range
        end
    end 
    delta_t(k) = mean(diff(datk.ftm)); %time [s] for each data point
end

N_detach = hist(Fdetach,force_bins);
N_f = sum(Nkf,1);
eps = (1/mean(delta_t)).*N_detach./N_f;
jlowN = find(N_f<N_f_min);
eps(jlowN) = NaN;

j_plus = find(force_bins>0);
j_minus = find(force_bins<0);
j_plus_fit = find(force_bins>0 & abs(force_bins)<=Fmax_fit);
j_minus_fit = find(force_bins<0 & abs(force_bins)<=Fmax_fit);

f_eps_plus = fit_single_exponential(force_bins(j_plus_fit),eps(j_plus_fit))
f_eps_minus = fit_single_exponential(-force_bins(j_minus_fit), eps(j_minus_fit)) 
jplot=find(eps ~=0);
figure, plot(force_bins(jplot),eps(jplot),'o','linewidth',2,'MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerSize',15)
hold on, plot(force_bins(j_plus_fit), feval(f_eps_plus, force_bins(j_plus_fit)), 'Color',color,'LineWidth',2)
hold on, plot(force_bins(j_minus_fit),feval(f_eps_minus, -force_bins(j_minus_fit)),'Color',color,'LineWidth',2)
xlabel('Force (pN)'), ylabel('Unbinding Rate (1/s)')
set(gca,'yscale','log')

figure, hold on
hp = bar(force_bins(j_plus),N_detach(j_plus)./sum(N_detach),1);
hm = bar(-force_bins(j_minus), -N_detach(j_minus)./sum(N_detach),1);
hp.FaceColor = color;
hp.FaceAlpha = alpha;
%hp.FaceColor = [0.65, 0.65, 0.65]; %gray
%hp.FaceAlpha = 1;
%hp.FaceColor = [0.47, 0.67, 0.19]; %green
%hp.FaceAlpha = 0.5;
hp.EdgeColor = 'None';
hm.FaceColor = hp.FaceColor;
hm.FaceAlpha = hp.FaceAlpha;
hm.EdgeColor = hp.EdgeColor;
%hm.FaceColor = [0.47, 0.67, 0.19];
%hm.FaceAlpha = 0.5;
%hm.EdgeColor = 'None';
xlabel('Force (pN)')
ylabel('Fraction of Events')


