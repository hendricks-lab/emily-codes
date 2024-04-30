%Emily Prowse's code for reading intensity and arborization lengths
clear all; close all; clc; 

addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Neuron_htt_int/';

set(0,'DefaultFigureWindowStyle','docked');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});

cmap = colormap(hot(100));
colour_30Q=[0.5 0.5 0.5];
colour_45Q=cmap(20,:);%[0.75 0 0];
colour_65Q=cmap(25,:);%[0.5 0 0];
colour_81Q=cmap(30,:);%[0.25 0 0];

htt_30Q="30Q";
htt_45Q="45Q";
htt_65Q="65Q";
htt_81Q="81Q";

a=0;
b=0;
c=0;
d=0;
e=0;
f=0;
g=0;
h=0;

directory='/Users/emilyprowse/Documents/McGill/Hendricks_lab/Analysis/';
addpath(directory);

cd(directory);
ratio_int_file='*intensity*HTT*.csv';
arb_len_file='*extension*HTT*.csv';
fname_int=dir(fullfile(directory,ratio_int_file));
fname_arb=dir(fullfile(directory,arb_len_file));

fname_int.name;
data=readcell(fname_int.name); %load everything in data
idx_30Q=find(data(:,3)=="30Q");
Q30_int=cell2mat(data(idx_30Q,9));
idx_45Q=find(data(:,3)=="45Q");
Q45_int=cell2mat(data(idx_45Q,9));
idx_65Q=find(data(:,3)=="65Q");
Q65_int=cell2mat(data(idx_65Q,9));
idx_81Q=find(data(:,3)=="81Q");
Q81_int=cell2mat(data(idx_81Q,9));

fname_arb.name;
arb_dat=readcell(fname_arb.name); %load everything in data
idx_30Q_arb=find(arb_dat(:,3)=="30Q");
Q30_arb=cell2mat(arb_dat(idx_30Q_arb,8));
idx_45Q_arb=find(arb_dat(:,3)=="45Q");
Q45_arb=cell2mat(arb_dat(idx_45Q_arb,8));
idx_65Q_arb=find(arb_dat(:,3)=="65Q");
Q65_arb=cell2mat(arb_dat(idx_65Q_arb,8));
idx_81Q_arb=find(arb_dat(:,3)=="81Q");
Q81_arb=cell2mat(arb_dat(idx_81Q_arb,8));

   
int{1}=Q30_int;
int{2}=Q45_int;
int{3}=Q65_int;
int{4}=Q81_int;

arb{1}=Q30_arb;
arb{2}=Q45_arb;
arb{3}=Q65_arb;
arb{4}=Q81_arb;

all_mean_int{1}=mean(int{1},"omitnan");
all_mean_int{2}=mean(int{2},"omitnan");
all_mean_int{3}=mean(int{3},"omitnan");
all_mean_int{4}=mean(int{4},"omitnan");

all_mean_arb{1}=mean(arb{1},"omitnan");
all_mean_arb{2}=mean(arb{2},"omitnan");
all_mean_arb{3}=mean(arb{3},"omitnan");
all_mean_arb{4}=mean(arb{4},"omitnan");

%Plotting the mean intensity for each cell
figure('Name','Mean_Intensity_Per_Cell','NumberTitle','off'), hold on,
plot(1,all_mean_int{1},'k_','MarkerSize',30);
plot(2,all_mean_int{2},'k_','MarkerSize',30);
plot(3,all_mean_int{3},'k_','MarkerSize',30);
plot(4,all_mean_int{4},'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(int{1}),1),int{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(int{2}),1),int{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(int{3}),1),int{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(int{4}),1),int{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
xlim([0 5]);
xticks([1 2 3 4]);
xticklabels({'30Q','45Q', '65Q', '81Q'})
ylabel('Htt/Tubb3 Intensity');
publication_fig(0,0,1);
box("off");
set(gca,'FontSize',40);

threshold=200;
arb_above_thres{1}=arb{1}(find(arb{1}>threshold));
arb_above_thres{2}=arb{2}(find(arb{2}>threshold));
arb_above_thres{3}=arb{3}(find(arb{3}>threshold));
arb_above_thres{4}=arb{4}(find(arb{4}>threshold));

arb_below_thres{1}=arb{1}(find(arb{1}<=threshold));
arb_below_thres{2}=arb{2}(find(arb{2}<=threshold));
arb_below_thres{3}=arb{3}(find(arb{3}<=threshold));
arb_below_thres{4}=arb{4}(find(arb{4}<=threshold));

figure('Name','Arborization_lengths_thresholded','NumberTitle','off'), hold on,
plot(1,all_mean_arb{1},'k_','MarkerSize',30);
plot(2,all_mean_arb{2},'k_','MarkerSize',30);
plot(3,all_mean_arb{3},'k_','MarkerSize',30);
plot(4,all_mean_arb{4},'k_','MarkerSize',30);
scatter(1-0.125+0.25*rand(numel(arb_below_thres{1}),1),arb_below_thres{1},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(arb_below_thres{2}),1),arb_below_thres{2},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(arb_below_thres{3}),1),arb_below_thres{3},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(arb_below_thres{4}),1),arb_below_thres{4},25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(1-0.125+0.25*rand(numel(arb_above_thres{1}),1),(threshold+50)*(ones(numel(arb_above_thres{1}),1)),25,'MarkerFaceColor','none','MarkerEdgeColor',colour_30Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(2-0.125+0.25*rand(numel(arb_above_thres{2}),1),(threshold+50)*(ones(numel(arb_above_thres{2}),1)),25,'MarkerFaceColor','none','MarkerEdgeColor',colour_45Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(3-0.125+0.25*rand(numel(arb_above_thres{3}),1),(threshold+50)*(ones(numel(arb_above_thres{3}),1)),25,'MarkerFaceColor','none','MarkerEdgeColor',colour_65Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);
scatter(4-0.125+0.25*rand(numel(arb_above_thres{4}),1),(threshold+50)*(ones(numel(arb_above_thres{4}),1)),25,'MarkerFaceColor','none','MarkerEdgeColor',colour_81Q,'MarkerEdgeAlpha',0.5,'LineWidth',2);

xlim([0 5]);
xticks([1 2 3 4]);
xticklabels({'30Q','45Q', '65Q', '81Q'})
ylabel('Arborization Length (\mum)');
publication_fig(0,0,1);
box("off")
set(gca,'FontSize',40);

%Statistical testing
[int_test, p_int_test]=ttest2(int{1},int{2});
[int_test_2, p_int_test_2]=ttest2(int{1},int{3});
[int_test_3, p_int_test_3]=ttest2(int{1},int{4});


[arb_test, p_arb_test]=ttest2(arb{1},arb{2});
[arb_test_2, p_arb_test_2]=ttest2(arb{1},arb{3});
[arb_test_3, p_arb_test_3]=ttest2(arb{1},arb{4});

titles={'30Q vs 45Q';'30Q vs 65Q';'30Q vs 81Q'};
% int_stats=[p_int_test;p_int_test_2;p_int_test_3];
% int_arb_stats=[p_int_arb_test;p_int_arb_test_2;p_int_arb_test_3];
%ANOVA
int_values=[int{1};int{2};int{3};int{4}];
int_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(int{1}) numel(int{2}) numel(int{3}) numel(int{4})])]';
int_data=table(int_values,int_groups);

[int_p,int_tbl_anova,int_stats]=anova1(int_values,int_groups,"on");

[int_results,~,~,gnames]=multcompare(int_stats);

int_tbl = array2table(int_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
int_tbl.("Group A") = gnames(int_tbl.("Group A"));
int_tbl.("Group B") = gnames(int_tbl.("Group B"));


arb_values=[arb{1};arb{2};arb{3};arb{4}];
arb_groups=[repelem([htt_30Q, htt_45Q, htt_65Q, htt_81Q], [numel(arb{1}) numel(arb{2}) numel(arb{3}) numel(arb{4})])]';
arb_data=table(arb_values,arb_groups);

[arb_p,arb_tbl_anova,arb_stats]=anova1(arb_values,arb_groups,"on");

[arb_results,~,~,gnames]=multcompare(arb_stats);

arb_tbl = array2table(arb_results,"VariableNames",["Group A","Group B","Lower Limit","Difference","Upper Limit","P-value"]);
arb_tbl.("Group A") = gnames(arb_tbl.("Group A"));
arb_tbl.("Group B") = gnames(arb_tbl.("Group B"));

fileprefix_stats='20231105_htt_tub3_arb_stats_';
% 
tempdir_1 = save_dir;   % Your destination folder
FolderName_1 = tempdir_1;   % Your destination folder
writetable(int_tbl,fullfile(FolderName_1, [fileprefix_stats,'int.csv']),'WriteRowNames',true);
writetable(int_data,fullfile(FolderName_1, [fileprefix_stats,'int_data.csv']),'WriteRowNames',true);
writetable(arb_tbl,fullfile(FolderName_1, [fileprefix_stats,'arb.csv']),'WriteRowNames',true);
writetable(arb_data,fullfile(FolderName_1, [fileprefix_stats,'arb_data.csv']),'WriteRowNames',true);

fileprefix='20231105_htt_tub3_arb_';

FolderName = save_dir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
end