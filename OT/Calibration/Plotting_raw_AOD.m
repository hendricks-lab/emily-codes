%testing TF input file
clear
%close all
clc
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/oras_codes-master/Optical Trap Calibration/')
%%% CHECK Fexc!!!%%%%
fileslocation = '/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/'; % defining the path to where the transferfunction.mat files are
%%% CHECK Fexc!!!%%%%
cd(fileslocation) %change the cd (current directory) to where the files are
location_in_numbers = double(fileslocation);
%dinfo = dir('*.txt'); %dir() shows all the files with .txt in the cd
dinfo = dir('20240416_kolf4_OT_s3_cal_5.txt'); %
Amp_wave=1;
beta=2000;%from calibration nm/V
cg = 1740;%1730;  %Ora 18/03/2021
fsamp_in = 20000;
% Columns in text file[QUADy input, QUADx input,SumSignal,AODx_input]
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
Fexc_in = [0.023  0.0357  0.053   0.097   0.17  0.27    0.43    0.57  0.813   1.23   2.11  3.19  5.07     8.77    19.9     36.7     83       163      349.1    547.1     899       1347];
for nfile = 1:length(dinfo)
thisfilename = dinfo(nfile).name; % currently selected file name in characters
    thisfilename = double(thisfilename); % currently selected file name in numerical values
    %%%%% HERE IS WHERE THE FUNCTION WOULD GO
    % - FUNCTION
    %Read file (Within the function)
    datdir = char(location_in_numbers); %converting numerical inputs back to characters
    filename = char(thisfilename); %converting numerical inputs back to characters
%     file name = strrep(filename,'.',''); % this replaces the "." of the original file name, to "" (blank). Weird glitch in the file name!!
    addpath(datdir);
    fileID = fopen(filename);
rawdata = textscan(fileID, '%f %f %f %f', 'HeaderLines', 0);
    xQUAD = rawdata{QUADx_input}.*kpd; % [V]
    yQUAD = rawdata{QUADy_input}.*kpd; % [V]
    sumQUAD_mean = mean(rawdata{3}).*kpd; % [V] % just to give an idea of the sum.
    xAOD = (rawdata{AODx_input}-75).*cg; % [nm]
    xAOD_v=rawdata{AODx_input}-75;
    Ldata = length(xAOD);
    Amp=Amp_wave*cg;
    time=[1:Ldata].*(1/fsamp_in);
% figure
% %plot(time,[xQUAD-mean(xQUAD)]*beta,'k'), hold on
% %plot(1:Ldata,yQUAD,'b'), hold on
% plot(time,xAOD,'r'), hold on
% legend(['Raw, Amp= ' num2str(Amp) 'nm'])
figure
plot(time,xQUAD-mean(xQUAD),'k'),hold on
plot(time,xQUAD,'r'),hold on
plot(time,xAOD_v,'b')
% legend('QPD V',['Raw, Amp= ' num2str(Amp_wave) 'V'])

figure
plot(time,yQUAD,'r'),hold on
end
