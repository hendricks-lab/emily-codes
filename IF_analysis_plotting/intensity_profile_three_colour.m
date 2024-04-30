%% Emily Prowse
%October 14, 2022
%Function for making 3 colour intensity line profiles from cilia images.
%Need to work on getting it to plot the intensities from all images. 

addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

clear all;close all; clc;

set(0,'DefaultFigureWindowStyle','docked');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});

save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Cilia_line_profiles/';

distance=[];
intensity=[];

for k_choose=1:3

    if k_choose==1 %Red
    cd('/Volumes/Emily_htt_2/Cilia/20221210/');
    IF_analysis_file='*linescan*r.csv';


    elseif k_choose == 2    % Green
    cd('/Volumes/Emily_htt_2/Cilia/20221210/');
    IF_analysis_file='*linescan*g.csv';

    
    elseif k_choose == 3   % Blue
    cd('/Volumes/Emily_htt_2/Cilia/20221210/');
    IF_analysis_file='*linescan*b.csv';


    end

fname=dir(IF_analysis_file);

for k=1:length(fname)
    
    data=readmatrix(fname(k).name); %load everything in the file
    distance{k,k_choose}=data(1:end, 1); 
    intensity=(data(1:end, 2));
    rel_intensity{k,k_choose}=intensity./(max(intensity));
    
end
end 

for m=1:length(distance)
    
    figure('Name',sprintf('tri_colour_intensity_profile %i', m),'NumberTitle','off'), hold on,
%     figure(m)
    plot(distance{m,1},rel_intensity{m,1},'-r'), hold on;
    plot(distance{m,2}, rel_intensity{m,2}, '-g'), hold on;
    plot(distance{m,3}, rel_intensity{m,3}, '-b');
    xlabel('Distance (\mum)');
    ylabel('Intensity');
    publication_fig(0,0,1);
end


% figure('Name','tri_colour_intensity_profile','NumberTitle','off'), hold on,


fileprefix='20221210_cilia_lineprofiles_';

FolderName = save_dir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end