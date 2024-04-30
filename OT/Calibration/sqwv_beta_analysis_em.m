%26 Mar 2021        mcgill      adam g. hendricks
% analyze square wave excitation data to estimate beta (nm/V)

% Commentary by Ora Cohen (ora.cohen@mail.mcgill.ca)

% Original code name: trap_sensitivity_step_excitation.m

%29.03.21 - Ora's edits
%           Changed fsamp to 2e4 vs 2e3 (bc it's actually 2e4)
%           Changed Cg from 1750 to 1730 bc I'm cocky like that
%           Added missing factor of 2 to uncorrected beta equation (Beta)
%           Added beta from base-to-spike signal (beta_base_avg)
%          (account for possibility of bead not reaching the trap center)

%05.04.2021 - More of Ora's edits
%             Corrected Beta2 to use dat0up and dat0dwn (instead of up twice) 
%             Corrected base-to-spike equation
%             Added Beta2 base-to-spike version

clear all, close all

Cg = 1730; %[nm/MHz] calibration constant for AOD (aka cx), named cg for curious george completely arbirarily by Ora and Adam collectively 

datdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/20210416/'; 
fname = '20210416_cell_beta_amp_70nm_6.txt'; %'20210408_mac_nd_beta_11.txt';
Fsamp = 2e4; %[Hz] sampling frequency
cht=0;
chxpd=1; %channel (text column) of photodiode x signal (on Labview, called y signal in other codes)
Gxpd=1;
chypd=2; %channel (text column) of photodiode y signal (on Labview, called x signal in other codes)
Gypd=1;
chxaod=4; %channel (text column) of aod x signal
Gaod=1;
chxstage=0;
Gxstage=1;

dat = textread(fullfile(datdir,fname),'','headerlines',1); %reading text file
dat0 = dat - mean(dat,1); %getting mean-centered data
time = (1/Fsamp).*[0:numel(dat(:,1))-1]'; %finding time in seconds

jp = find(time>0 & time<10); %used to plot signal over time
figure(1), plot(time(jp),dat0(jp,chypd),time(jp),dat0(jp,chxaod)) 

%divide signal into right and left steps
jup = find(dat0(:,chxaod)>0); %gives index of all data where aod is positive
jdwn = find(dat0(:,chxaod)<0); %gives index of all data where aod is negative

%find starts of each step
jz_up = find(diff(jup)>1); %gives index of jup where difference between two values is > 1 (a sudden step)
jz_dwn = find(diff(jdwn)>1); %gives index of jdwn where difference between two values is > 1 (a sudden step)

% if jup(1)==1,
%     jz_up=jz_up(2:end);
% end
% if jdwn(1)==1,
%     jz_dwn=jz_dwn(2:end);
% end

% average all steps

%find how long a single step is (num data points between each step)
Lup = jup(jz_up(2)+1)-jup(jz_up(2)); %length of a single step (index of the start of the next step - index of last step)
dat0up=zeros(Lup,4); %setting up a matrix to hold averaged data over all "up" steps
Nzup=numel(jz_up)-1; %number of "up" steps
dat0dwn = zeros(Lup,4); %setting up a matrix to hold averaged data over all "down" steps
Nzdwn = numel(jz_dwn)-1; %number of "down" steps

%
%averaging the data over all "up" steps
for kz = 1:Nzup,
    dat0up=dat0up+dat0(jup(jz_up(kz)):jup(jz_up(kz)+1)-1,:);
end
dat0up=dat0up./Nzup;

%averaging the data over all "down" steps
for kz = 1:Nzdwn,
    dat0dwn = dat0dwn + dat0(jdwn(jz_dwn(kz)):jdwn(jz_dwn(kz)+1)-1,:);
end
dat0dwn = dat0dwn./Nzdwn;

%plot averaged data from "up" and "down" steps
figure, plot(time(1:Lup),dat0up(:,chypd))
hold on, plot(time(1:Lup),dat0dwn(:,chypd))

%only fit exponential to initial signal (cut off first few points bc
%signal was still increasing from slight delay in qpd response to trap mvmt)
jfit = find(time>0.0002 & time < 0.01);

%fit exponential to the averaged qpd "up" signal
fup = fit(time(jfit), dat0up(jfit,chypd), 'exp2');
hold on, plot(time(jfit),feval(fup,time(jfit)),'linewidth',2)
V_0up=fup.a+fup.c; %QPD voltage immediately following an "up" step
A_0up=mean(dat0up(2:end,chxaod)); %AOD signal (amplitude)
%phat = mle(dat0up(10:end,chypd)+,'distribution','exp')

%fit exponential to the averaged qpd "down" signal
fdwn = fit(time(jfit), dat0dwn(jfit,chypd), 'exp2');
hold on, plot(time(jfit),1*feval(fdwn,time(jfit)),'linewidth',2)
V_0dwn = fdwn.a + fdwn.c; %QPD voltage immediately following a "down" step
A_0dwn = mean(dat0dwn(2:end,chxaod)); %AOD signal (amplitude)

%baseline QPD signal (after bead has relaxed following a step)
V_basedwn = (mean(dat0dwn(end-floor(Lup./10):end,chypd)));
V_baseup = (mean(dat0up(end-floor(Lup./10):end,chypd)));

%this beta is based on a kpd value of 1, not 10! Also does not account
%for bead not reaching the center of the trap between steps
disp('Based on exp fit spike-to-spike (kpd = 1)')
Beta = 2*Cg*(A_0dwn-A_0up)/(V_0up-V_0dwn) %Beta based on exponential fit (with corrected factor of 2)
    
Beta_baseup = abs(Cg*(A_0dwn-A_0up)/(V_0up-V_basedwn)); %Beta correcting for initial signal =! 0
Beta_basedown = abs(Cg*(A_0dwn-A_0up)/(V_0dwn-V_baseup)); %Beta correcting for initial signal =! 0

%This is a corrected beta to account for the bead not reaching the trap
%center between steps. Based on a kpd of 1
disp('Based on exp fit base-to-spike (kpd = 1)')
beta_base_avg = (Beta_baseup+Beta_basedown)./2 %Beta based on exponential fit corrected for non-zero base signal

disp('Based on signal max spike-to-spike (kpd = 1)')
Beta2 = 2.*Cg.*(A_0dwn-A_0up)/(max(abs(dat0up(:,chypd)))+max(abs(dat0dwn(:,chypd)))) %Beta based on signal max value only, not exponential fit (with corrected factor of 2)

Beta2Up = abs(2*Cg*A_0up/((max(dat0up(:,chypd))-V_basedwn)));
Beta2Down = abs(2*Cg*A_0dwn/((min(dat0dwn(:,chypd))-V_baseup)));

disp('Based on signal max base-to-spike (kpd = 1)')
Beta2Ora = (Beta2Up+Beta2Down)/2 %Beta based on signal max corrected for non-zero base signal

%plot calibrated 
figure, plot(time,Cg.*dat0(:,chxaod),'r', 'linewidth', 2)
hold on, plot(time, Beta.*dat0(:,chypd)+Cg.*dat0(:,chxaod),'b', 'linewidth', 2)
title('Calibrated bead position based on uncorrected beta')
xlabel('Time (s)'), ylabel('position (nm)'), legend('trap','bead')

%plot calibrated 
figure, plot(time,Cg.*dat0(:,chxaod),'r', 'linewidth', 2)
hold on, plot(time, beta_base_avg.*dat0(:,chypd)+Cg.*dat0(:,chxaod),'b', 'linewidth', 2)
title('Calibrated bead position based on base-corrected beta')
xlabel('Time (s)'), ylabel('position (nm)'), legend('trap','bead')
