%% Emily Prowse
%October 13, 2022
%Function for making box plots from csv files with information on cell
%morphology in this example

addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

clear all;close all; clc;

set(0,'DefaultFigureWindowStyle','docked');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});

save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Morphology/';

for k_choose=1:2

    if k_choose == 1    % WT 
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Morphology/');
    per_cell_num=[];
    per_cell_area=[];
    IF_analysis_file='*WT*morphology.csv';

    elseif k_choose == 2    % S421D
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Morphology/');
    per_cell_num=[];
    per_cell_area=[];
    IF_analysis_file='*S421D*morphology.csv';

    end
    
fname=dir(IF_analysis_file);
data=readmatrix(fname.name); %load everything in the file
 
area{k_choose}=data(1:end,2); % Area
mean_int{k_choose}=data(1:end,3); % Mean actin intensity
circ{k_choose}=data(1:end,4); % Circularity
ar{k_choose}=data(1:end,5); % Aspect Ratio
rnd{k_choose}=data(1:end,6); % Roundness
sldty{k_choose}=data(1:end,7); % Solidity
    
end

%Plot the data as histograms
figure('Name','Cell_area_hist','NumberTitle','off'), hold on,
histogram((area{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((area{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Cell area (\mum^2)');
ylabel('Frequency');
publication_fig(0,0,1);

figure('Name','Actin_intensity_hist','NumberTitle','off'), hold on,
histogram((mean_int{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((mean_int{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Actin Intensity');
ylabel('Frequency');
publication_fig(0,0,1);

figure('Name','Circularity_hist','NumberTitle','off'), hold on,
histogram((circ{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((circ{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Circularity');
ylabel('Frequency');
publication_fig(0,0,1);

figure('Name','Aspect_ratio_hist','NumberTitle','off'), hold on,
histogram((ar{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((ar{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Aspect Ratio');
ylabel('Frequency');
publication_fig(0,0,1);

figure('Name','Roundness_hist','NumberTitle','off'), hold on,
histogram((rnd{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((rnd{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Roundness');
ylabel('Frequency');
publication_fig(0,0,1);

figure('Name','Solidity_hist','NumberTitle','off'), hold on,
histogram((sldty{1}),20,'facecolor',[0.5 0.5 0.5], 'edgecolor','none','facealpha',0.5), hold on,
histogram((sldty{2}),20, 'facecolor',[0 0.4 0.4], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Solidity');
ylabel('Frequency');
publication_fig(0,0,1);


figure('Name','Cell_area_violin','NumberTitle','off'), hold on,
violin(area, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(area{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(area{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Cell area (\mum^2)');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Actin_intensity_violin','NumberTitle','off'), hold on,
violin(mean_int, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(mean_int{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(mean_int{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Actin Intensity');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Circularity_violin','NumberTitle','off'), hold on,
violin(circ, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(circ{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(circ{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on'); 
ylabel('Circularity');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Aspect_ratio_violin','NumberTitle','off'), hold on,
violin(ar, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(ar{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(ar{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Aspect Ratio');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Roundness_violin','NumberTitle','off'), hold on,
violin(rnd, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(rnd{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(rnd{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Roundness');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Solidity_violin','NumberTitle','off'), hold on,
violin(sldty, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5;0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(sldty{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(sldty{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
publication_fig(0,0,1);
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',40);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Solidity');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

%Statistical test suggested by Adam
[area_test, p_area_test]=ttest2(area{1},area{2});
[int_test, p_int_test]=ttest2(mean_int{1},mean_int{2});
[circ_test, p_circ_test]=ttest2(circ{1},circ{2});
[ar_test, p_ar_test]=ttest2(ar{1},ar{2});
[rnd_test, p_rnd_test]=ttest2(rnd{1},rnd{2});
[sldty_test, p_sldty_test]=ttest2(sldty{1},sldty{2});

% fileprefix='20221013_WT_vs_S421D_morphology_';
% 
% FolderName = save_dir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end