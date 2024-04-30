clear all %, close all, warning off
close all %Emily added
clc;
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/loic_codes-master 2/adams-method/') %path where some functions are stored

datdir = '/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/';
% datdir = '/Volumes/LACIE-500GB/emily-backup/20220122_30Q_DIV7_OT';
% datdir = '/Volumes/Emily_htt_2/HEK_OT/20210831_S421D_OT_early_late/'
caldir='/Volumes/Emily_2022/Neurons/Good_OT_Calibrations/';
addpath(datdir)

%set filename to multiharmonic measurement
flist{1} = {'18.04.24-11.53-TF-20240416_kolf4_OT_s3_cal_6.mat'}; %name of file to load
filename =  fullfile(datdir,flist{1}{1}); % selecting the file, using the file path and file name
data=load(filename);

    cd(caldir)

    flist=dir([caldir,'*.mat']);

    Temp = 37; % [deg C]
    kT = 1.381e-2*(Temp+273.15); % [K]
    Rbead = 100; % Radius of bead [nm], Emily modified for nanodiamonds (not sure if this is even being used)
    Zbead = 2000; % Height of bead [nm] (not currently being used (August 2019) but it's something that could be incorporated.
    falias=20000;
    W_vector = ones(1,length(data.Fexc)); 
    H=data.H;
    PHfrf=data.PHfrf;
    Wmag = 5e-2; %Weight applied to magnitude (from normal code)
    Wphase = 0.05.*Wmag; %Weight applied to phase (from normal code)
    W = Wmag.*W_vector;
    Wphase = Wphase.*W_vector;
    fsamp=data.fsamp;

for n=1:numel(flist)
    u0_loaded=load(flist(n).name);
    close()
    close()
    close()

    fdata1=u0_loaded.fdata1;
    ydat10=u0_loaded.ydat10;

    xdata = [data.Fexc,data.Fexc,fdata1'];
%     xdata2 = xdata(find(xdata>100 & xdata <2000));
    ydata = [zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
%     ydata2 = ydata(find(xdata>100 & xdata <2000));
    
    PYm=2*mean(diff(data.FP_y)).*data.PYm_y;

    if data.C(:,14:end) > 0.8
    
%     figure(n*100),
%     plot(data.FP_y, PYm,'b-'), hold on
%     plot(u0_loaded.freqplot, u0_loaded.PS_fit, 'c-')
%     set(gca, 'YScale', 'log');
%     set(gca, 'XScale', 'log');
%     xlabel('Frequency (Hz)');
%     ylabel('Power (V^2s');

    

    residual_indices=vertcat(find(data.FP_y>=20 & data.FP_y<=80), find(data.FP_y>=300 & data.FP_y<=1500));

    figure(n*100),
    plot(data.FP_y, PYm,'b-'), hold on
    plot(u0_loaded.FP_y, u0_loaded.PYm, 'm-')
    plot(data.FP_y(residual_indices),PYm(residual_indices),'y-');
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('Power (V^2s');

    residual_raw_y=PYm(residual_indices)-u0_loaded.PYm(residual_indices);
    residual_raw_x=data.FP_y(residual_indices);

    mean_res_raw{n}=mean(residual_raw_y);

%     figure(n*100000)
%     plot(residual_raw_x,residual_raw_y,'m-');
    
    lb=[0.8*u0_loaded.beta     0.1*u0_loaded.gammar    0.95*u0_loaded.alpha     0.25.*u0_loaded.k_trap   0.1*u0_loaded.k_cyt0    0.1*u0_loaded.k_cyt1    1e-1    0]; %Emily
    ub=[1.5*u0_loaded.beta     10*u0_loaded.gammar     1.05*u0_loaded.alpha     5*u0_loaded.k_trap     10.*u0_loaded.k_cyt0    10*u0_loaded.k_cyt1    1e1     1.01]; %Emily
    
    u0=[u0_loaded.beta, u0_loaded.gammar,  u0_loaded.alpha, u0_loaded.k_trap, u0_loaded.k_cyt0, u0_loaded.k_cyt1, 10, 1];

    options_fit1 = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-3);
    [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_r9,u0,xdata,ydata,lb,ub,options_fit1,abs(H(:))',PHfrf,W,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias); 

    residi{n}=residual;
    resnorm_all{n}=resnorm;
%     errfit = theor_trap_frf_r9(u2,xdata,abs(H(:))',PHfrf,W2,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias); % this function gives the weighted residual
    else
        disp('Coherence too low to fit')
    end

  
end


% for k=1:numel(residi)
% %     mean_res{k}=mean(residi{k});
% %     figure(k*1000),
% %     plot(residi{k}(1:200));
%     mean_res_raw{k}=mean(cell2mat(residual_raw_y{k}));
% %     figure(10000)
% %     plot(cell2mat(residual_raw_x{k}), cell2mat(residual_raw_y{k}))
% end

% minimum=min(abs(cell2mat(mean_res)));
minimum=min(abs(cell2mat(mean_res_raw)));



% best_cal_idx=find((abs(cell2mat(mean_res))==minimum));
best_cal_idx=find((abs(cell2mat(mean_res_raw))==minimum));
    