%%% Loïc Chaubet
%%% March 2019
%%% 
%%% 
%%% This code runs:
%%% function [H,Fexc,C,FP,PYm,Pyy,Pyy2,fsamp,length_used_out,sumQUAD_mean,Npeaks_out,all_loc_out,K_int_out,frequency_range_tolerance]= Step1_TF_Loic_2018b_FUNCTION(N_channels,kpd,AODx_input,QUADx_input,nfftcoef,truncate_to_period,hannw,plot_range,cut_peaks_data,mean_threshold,cut_threshold,mov_wind_s,plot_cut,Fexc_in,fsamp_in,cg,fileID)
%%% on ALL .txt files in the folder.

%%% Original name: Step1_TF_Loic_2018b_read_all_using_function

% Ora's changes 
%     Set constants to correct values instead of asking user for input (10/2020)
%     Added excitation frequencies from new wave (12/2020)
%     Added y power spectrum calculation/output (03/2021)
%     Changed Cg to 1730 (03/2021)
%     cleaned up some extra lines (04/2021)
%     Changed kpd to 1 (04/2021)


clear
%close all
clc
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/loic_codes-master 2/TF/')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/oras_codes-master/Optical Trap Calibration/')

%%% CHECK Fexc!!!%%%%

fileslocation = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/'; % defining the path to where the transferfunction.mat files are
% fileslocation = '/Volumes/Emily_htt_2/20210810_OT_WT_fn_early_to_late/'; % defining the path to where the transferfunction.mat files are
%%% CHECK Fexc!!!%%%%

cd(fileslocation) %change the cd (current directory) to where the files are
% location_in_numbers = double(fileslocation);

%dinfo = dir('*.txt'); %dir() shows all the files with .txt in the cd

% Emily's code to make fake data for OT to test out the frequency response
t=linspace(1,275,5500000);
input_frequency=4*sin(2*pi*t);
output_frequency_in_phase=2*sin(2*pi*t);
output_frequency_out_of_phase=2*cos(2*pi*t);

% plot(t,input_frequency), hold on
% plot(t,output_frequency_in_phase), hold on
% plot(t,output_frequency_out_of_phase);

dinfo = dir('*fake*.txt'); %Emily modified for fake data

kpd = 1; %input('What is the voltage constant multiplier? (10 for recent data, 1 for older data) '); %except apparently now it's 1 again?
AODx_input = 4; %input('In which channel is the AODx input signal? (4 for recent data, 5 for older data) ');
QUADx_input = 2; %input('In which channel is the QUADx measured signal? (2 for recent data, ?? for older data) ');
QUADy_input = 1; %input('In which channel is the QUADy measured signal? (1 for recent data, ?? for older data) ');
N_channels = 4; %input('How many channels total? (4 for recent data, 6 for older data) ');
nfftcoef = 1.0; %input('What is the coefficient of NFFT? (1.0 typically, +/- 10% max...need to keep it consistant though!!) ');
cutpeaks = 1; % input('Any data to be cut, due to joystick movement or major drift events? (yes = 1, no = 0) ');
manualcut = 0; % cut data?
manualcut_data_end = 0; % cut last 20 seconds of data

plot_cut = 0; % plot cut?
plot_range = 0; % plot range (TF estimate's magnitude and coherence near the Fexc's)
hannw = 0; % hannw = 1 will use hann window. Set it to 0 to use Hamming window

cg = 1730; %Ora 18/03/2021
fsamp_in = 20000;
%range_tolerance = 1

% 29 November 2017, Dec 12 to end of February (1.0x)
      %  Fexc_in =         [0.023   0.053   0.17   0.43   1.37    2.37    5      7     19     27     67     87     157    337      899      1347      2441];

% March 1st 2018 (1.0x frequencies) to May 11th 2018
%Fexc_in =  [0.023   0.053   0.097  0.17   0.43   1.37    2.37    5      9.1     19     37     87     157    337      899      1347      2441];

% May 11th to June 19th. NOTE: it seems after 21 freqs, LABVIEW cannot generate proper input signal!?, and there is crosstalk between 0.79 and 2.37!!!
%Fexc_in =  [0.023   0.053   0.097  0.17  0.27   0.43  0.79  1.37    2.37    5      9.1     19     37     87     157      337      547    899    1347];

% June 19th to july 10th
%Fexc_in =    [0.023   0.053   0.097  0.17  0.27   0.43  0.79  1.37    2.41    5.1      9.1     19     37     87     157      337      547    899    1347];

% july 10th to August 24th
%Fexc_in =              [0.023   0.053   0.097   0.17   0.43   0.57    1.23    2.11    5.07     8.77    19.9     36.7      83       163      349.1      899       1347];

% anything past August 24th 2018 using default wave - 4 columns of data only
% Fexc_in = [0.023  0.053   0.097   0.17  0.27    0.43    0.57     1.23   2.11    5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

% anything past August 24th 2018 using default wave - 4 columns of data only
Fexc_in = [1]; %Emily modified for fake data

% any measurements 1 minute long - 4 columns of data only, take out
% longest frequency
%Fexc_in = [0.053   0.097   0.17  0.27    0.43    0.57     1.23   2.11    5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

%  Ora's super fun 16/12/2020 experiment, 2 new freqs - 4 columns of data only
%Fexc_in = [0.023  0.037  0.053   0.097   0.17  0.27    0.43    0.57  0.91   1.23   2.11    5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

%  Ora's super fun 16/12/2020 experiment, 3 new freqs - 4 columns of data only
%Fexc_in = [0.023  0.037  0.053   0.097   0.17  0.27    0.43    0.57  0.91   1.23   2.11  3.03  5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

%  Ora's super fun 12/01/2021 experiment, 3 new freqs - 4 columns of data only
%Fexc_in = [0.023  0.0357  0.053   0.097   0.17  0.27    0.43    0.57  0.813   1.23   2.11  3.19  5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

% FOR TESTING
%Fexc_in = [0.023  0.053   0.097   0.17  0.27    0.43    0.57     1.23   2.11 ]%   5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];
%Fexc_in = [1347]%  0.053   0.097 0.17 0.27]% 
%Fexc_in = [19.9 ]
%Fexc_in = [0.17 0.27]%0.17  0.27    0.43    0.57  1.23 ];%     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];
%Fexc_in = [0.57 1.23  2.11    5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];

%NOTE: starting at 349.1 Hz, the bin Hz are slightly misaligned [349.1004  547.1001  899.0007  1347.0004]

if cutpeaks == 1
    std_threshold = 10;%input('Threshold for the ratio of standard deviation (typically 10) ');
    cut_threshold = 0.5;%input('Threshold for the size of peak to be cut (typicall 0.5) ');
    mov_wind_s  = 0.5;%input('Moving window length in seconds over which to integrate the full joystick movement/major drift (typically 0.5) ');
else
    std_threshold = [];
    cut_threshold = [];
    mov_wind_s = [];
end

if kpd ~= 1 && kpd ~=10
    disp('Script stopped: kpd must be 1 or 10')
    return
end

for nfile = 1:length(dinfo)
    
    thisfilename = dinfo(nfile).name; % currently selected file name in characters
    thisfilename = double(thisfilename); % currently selected file name in numerical values
  
    %%%%% HERE IS WHERE THE FUNCTION WOULD GO
    % - FUNCTION   
    %Read file (Within the function)
%     datdir = char(location_in_numbers); %converting numerical inputs back to characters
    filename = char(thisfilename); %converting numerical inputs back to characters
    filename = strrep(filename,'._',''); % this replaces the "._" of the original file name, to "" (blank). Weird glitch in the file name!!
%     addpath(datdir);
    fileID = fopen(filename);
    
    [H,Fexc,C,FP_x,PYm_x,FP_y,PYm_y,Pyy,Pyy2,fsamp,length_used,sumQUAD_mean,Npeaks,all_loc,K_integer,freq_tolerance,NFFT,number_windows,data_lost,Noverlap] = Step1_get_TF_PSx_PSy_FUNCTION_fakedata(N_channels,kpd,AODx_input,QUADx_input,QUADy_input,nfftcoef,hannw,plot_range,cutpeaks,std_threshold,cut_threshold,mov_wind_s,plot_cut,Fexc_in,fsamp_in,cg,fileID,manualcut,manualcut_data_end);
    
    Mfrf = 20.*log10(abs(H));
    PHfrf = unwrap(angle(H));
    PHfrf = abs(PHfrf - round(PHfrf./pi).*pi)*180/pi;% why??! 
    
    % plots the FRF and saves the handle to "hfrf"
    figure, subplot(211), semilogx(Fexc,Mfrf,'o','linewidth',2)
    ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
    subplot(212), semilogx(Fexc,PHfrf,'o','linewidth',2)
    ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
    hfrf=gcf;

    % plots the coherence and saves the handle to "hch"
    figure, semilogx(Fexc,C,'o','linewidth',2)
    ylabel('Coherence'), xlabel('f, Hz')
    hch=gcf;

    % plots the power spectrum and saves the handle to "hpsx"
    figure, loglog(FP_x,PYm_x,'b','linewidth',2), hpsx=gcf;
    xlabel('f (Hz)'), ylabel('Power Spectrum in X (V^2 s)')
    
    % plots the power spectrum and saves the handle to "hpsy"
    figure, loglog(FP_y,PYm_y,'b','linewidth',2), hpsy=gcf;
    xlabel('f (Hz)'), ylabel('Power Spectrum in Y (V^2 s)')
    
    % plots the different Welch's method parameters NFFT,number_windows,data_lost
    figure, subplot(311), semilogx(Fexc,NFFT,'o','linewidth',2),ylabel('NFFT')
    subplot(312), semilogx(Fexc, number_windows,'o','linewidth',2),ylabel('Number of windows')
    subplot(313), semilogx(Fexc, data_lost,'o','linewidth',2),ylabel('% of data discarded by Pwelch')
    hWelch=gcf;
    
    currentFolder = pwd; % registering the current folder name
    cd(currentFolder) % this makes sure the new .mat file is written in the current folder
   %
    if cutpeaks == 1
        %savename = sprintf('TF-MATLAB2018b-CUT-FUNCTION-%s',filename); % puts a bunch of things right before the file name string
%         figure(hfrf), title(['#peaks = ' num2str(Npeaks) ' at ' num2str((all_loc./fsamp)') ''])
%         figure(hch), title(['Cutting parameters: cut threshold = ' num2str(std_threshold) ' cut width threshold = ' num2str(cut_threshold) ' moving window size = ' num2str(mov_wind_s) ''])
    end
    %else
        %savename = sprintf('TF-%s',filename); % puts "TF-" right before the file name string
    %end
    %
    savename = sprintf('TF-%s',filename); % puts "TF-" right before the file name string
    
    date_now = datestr(now, 'dd.mm.yy-HH.MM');
    savename = [date_now,'-',savename]; % this adds the date and time at the end of the filename
    savename = strrep(savename,'.txt','.mat'); % this replaces the ".txt" of the original file name by ".mat"
    K_integer = K_integer;
    freq_tolerance = freq_tolerance;
    Noverlap;
    NFFT;

    save(savename) % this saves the whole workspace, including the 4 figures that were assigned to the 3 corresponding handles 
    close all
end

return

% raw = rawdata{QUADx_input}.*kpd;
% figure,
% plot(n,raw,'k'), hold on
% plot(n(indices_to_cut),raw(indices_to_cut),'g')
% legend('Raw')
% load(filename,'Fexc','H','FP','PYm','Pyy','Pyy2','fsamp','Rbead','Zbead','C','kT','Pxx');
