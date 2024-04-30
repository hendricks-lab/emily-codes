%Emily Prowse (updated Oct 29, 2021)

%This function is meant to measure fluorophore emission crosstalk
%(bleedthrough) or colocalization of two image channels.
%The images it takes in should be taken in two channels for a single
%fluorophore, the channel it is expected to emit in, and the channel that
%it was not intended to emit in but its fluorescence is observed. It will
%form a plot of the intensities registered for the two images, which should
%show a linear relation so that crosstalk can be measured and corrected for
%in ImageJ. 

%motor_colocalization
% To use this code, the directory window will appear when you click "run".
% You will need to select 3 things, the first is the directory where your
% images are stored, this is narrowed down by the pathname. Next, you will
% select the green image (a tiff file) , then the red image (a tiff file). 
%Note that the images must be grayscale for this code to run, if they 
%aren't use the rgb2gray() function.

addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

clear all;
close all;
clc;
% PathName = uigetdir('/Volumes/Emily/20210919/20210919_Tiffs/20210919_focus_zplane/') 
mdl_prev=[];
r_sqr_prev=[];
model=[];
r_squared=[];

% for condition = 1:4 %for khc/dynein
  for condition = 1:6 %for other kifs

mdl_prev=[];
r_sqr_prev=[];
% model=[];
% r_squared=[];
    
if condition == 1    % WT kif1a
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif1a/mt_reg/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif1a/kif1a/';    
elseif condition == 2    % S421D kif1a
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif1a/mt_reg/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif1a/kif1a/';
elseif condition == 3   % WT kif3a
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif3a/mt/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif3a/kif3a/'; 
elseif condition == 4   % S421D kif3a
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif3a/mt_reg/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif3a/kif3a/';
elseif condition == 5   % WT kinf16b
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif16b/mt_reg/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/WT/kif16b/kif16b/';
elseif condition == 6   % S421D kif16b
    folder_mt='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif16b/mt_reg/';
    folder_motor='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK293T_IF/S421D/kif16b/kif16b/';
end



fname_motor=dir(folder_motor);
fname_mt=dir(folder_mt);

thresholdedredintensity_previous=[];
thresholdedgreenintensity_previous=[];
dirsize=[];
% for k=1:length(fname_mt)
min(numel(fname_motor),numel(fname_mt));
if numel(fname_motor)==numel(fname_mt)
    dirsize=numel(fname_mt);
else
    display('Error: Directories are not the same size, using the directory with less files but you might want to check why this number is different')
    dirsize=min(numel(fname_motor), numel(fname_mt));
end

for k=1:dirsize
save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/IF_motor/';
close all


if startsWith(fname_motor(k).name,'.') || startsWith(fname_mt(k).name,'.')
    continue
end
                cd(folder_motor)
                Ig=imread(fname_motor(k).name); %read the green image data
                greenintensity = double(Ig(1:262144)); %sets y values as the green image intensity values
                %the number 262144 comes from the number of pixels in a
                %16-bit image, adjust based on your image size.
                yaxis = 'Microtubule Image Intensity';%sets label for the x axis (red image)
                cd(folder_mt)
%                 if startsWith(fname_mt(k).name,'.') || startsWith(fname_motor(k).name,'.')
%                 continue
%                 end
                Ir=imread(fname_mt(k).name); %read the red image data
                redintensity = double(Ir(1:262144)); %sets x values as the red image intensity values6
                xaxis = 'Motor Image Intensity'; %sets label for the y axis (green image)           
                
               %If a threshold is needed, here is the code to add it
                redintensityabovethreshold= find(redintensity>3000);
                greenintensityabovethreshold= find(greenintensity>200);
                
                indexforabovethreshold=redintensityabovethreshold.*ismember(redintensityabovethreshold, greenintensityabovethreshold);
                indexforabovethreshold(indexforabovethreshold==0)=[];
                thresholdedredintensity=redintensity(indexforabovethreshold);
                thresholded_red_all{condition}=thresholdedredintensity;thresholdedredintensity_previous;
                thresholdedgreenintensity=greenintensity(indexforabovethreshold);
                thresholded_green_all{condition}=thresholdedgreenintensity;thresholdedgreenintensity_previous;

%               %No threshold is required
%                 thresholdedredintensity=redintensity;
%                 thresholdedgreenintensity=greenintensity;                

                mdl = polyfit(thresholdedgreenintensity, thresholdedredintensity, 1); %create linear regression
                
                mdly= polyval(mdl,thresholdedgreenintensity); 
                %this calculates the y values using the polynomial fit model
                save_model=[mdl;mdl_prev];
                model{condition}=save_model;
                mdl_prev=save_model;
                
                correlation_coeff = corr2(thresholdedgreenintensity, thresholdedredintensity);
                r_sqr = power(correlation_coeff,2);
                save_rsqr=[r_sqr;r_sqr_prev];
                r_sqr_prev=save_rsqr;
                r_squared{condition}=save_rsqr;
                %this will print the value stored in model of the
                %polynomial fit, first the slope then the intercept value.
                
%                 figure('Name','Fraction of Microtubule Associated Kinesin Raw','NumberTitle','off'), hold on,
%                 plot(thresholdedgreenintensity, thresholdedredintensity, 'o', 'MarkerSize', 5, 'color', 'k'); %plotting the raw intensity values
%                 hold on %saves the previous plot and allows you to add the next plot to the same one.
%                 plot(thresholdedgreenintensity, mdly, 'color', 'b'); % plot the linear regression
%                 xlabel(xaxis); % x axis label, defined by variable xaxis
%                 ylabel(yaxis); % y axis label, defined by variable yaxis
% 
%                 figure('Name','Fraction of Microtubule Associated Kinesin Colormap','NumberTitle','off'), hold on,
%                 colormap(parula);
%                 Nbins=100;
%                 hold on, hb = binscatter(thresholdedgreenintensity, thresholdedredintensity);
%                 hb.NumBins = [Nbins Nbins];
%                 hb.ShowEmptyBins = 'off';
%                 hb.FaceAlpha = 0.75;
%                 xlabel(xaxis); ylabel(yaxis);
%                 colorbar;
%                 hold on
%                 plot(thresholdedgreenintensity, mdly, 'color', 'b', 'LineWidth', 1); % plot the linear regression
% 
%                 figure('Name','Fraction of Microtubule Associated Kinesin Colormap without linear regression','NumberTitle','off'), hold on,
%                 colormap(parula);
%                 Nbins=100;
%                 hold on, hb = binscatter(thresholdedgreenintensity, thresholdedredintensity); 
%                 hb.ShowEmptyBins = 'off';
%                 hb.FaceAlpha = 0.75;
%                 xlabel(xaxis); ylabel(yaxis);
%                 colorbar;
                
%                 fileprefix=append(fname_motor(k).name,'_thresholded_removedbad');
%                 FolderName = save_dir;   % Your destination folder
%                 FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
%                 for iFig = 1:length(FigList)
%                 FigHandle = FigList(iFig);
%                 FigName   = get(FigHandle, 'Name');
%                 savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%                 saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
                thresholdedredintensity=thresholdedredintensity_previous;
                thresholdedgreenintensity=thresholdedgreenintensity_previous;

                
% end
end                
  end

%Emily adding a portion to plot all the data from all the files of each condition together and fit it
% mdl_all_WT_khc = polyfit(thresholded_green_all{1}, thresholded_red_all{1}, 1); %create linear regression
% mdly_all_WT_khc= polyval(mdl_all_WT_khc,thresholded_green_all{1}); 
% %this calculates the y values using the polynomial fit model
% correlation_coeff_WT_khc = corr2(thresholded_green_all{1}, thresholded_red_all{1});
% r_sqr_WT_khc = power(correlation_coeff_WT_khc,2);
% 
% mdl_all_S421D_khc = polyfit(thresholded_green_all{2}, thresholded_red_all{2}, 1); %create linear regression
% mdly_all_S421D_khc= polyval(mdl_all_S421D_khc,thresholded_green_all{2}); 
% %this calculates the y values using the polynomial fit model
% correlation_coeff_S421D_khc = corr2(thresholded_green_all{2}, thresholded_red_all{2});
% r_sqr_S421D_khc = power(correlation_coeff_S421D_khc,2);
% 
% mdl_all_WT_dyn = polyfit(thresholded_green_all{3}, thresholded_red_all{3}, 1); %create linear regression
% mdly_all_WT_dyn= polyval(mdl_all_WT_dyn,thresholded_green_all{3}); 
% %this calculates the y values using the polynomial fit model
% correlation_coeff_WT_dyn = corr2(thresholded_green_all{3}, thresholded_red_all{3});
% r_sqr_WT_dyn = power(correlation_coeff_WT_dyn,2);
% 
% mdl_all_S421D_dyn = polyfit(thresholded_green_all{4}, thresholded_red_all{4}, 1); %create linear regression
% mdly_all_S421D_dyn= polyval(mdl_all_S421D_dyn,thresholded_green_all{4}); 
% %this calculates the y values using the polynomial fit model
% correlation_coeff_S421D_dyn = corr2(thresholded_green_all{4}, thresholded_red_all{4});
% r_sqr_S421D_dyn = power(correlation_coeff_S421D_dyn,2);
% %this will print the value stored in model of the
% %polynomial fit, first the slope then the intercept value.

mdl_all_WT_kif1a = polyfit(thresholded_green_all{1}, thresholded_red_all{1}, 1); %create linear regression
mdly_all_WT_kif1a= polyval(mdl_all_WT_kif1a,thresholded_green_all{1}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_WT_kif1a = corr2(thresholded_green_all{1}, thresholded_red_all{1});
r_sqr_WT_kif1a = power(correlation_coeff_WT_kif1a,2);

mdl_all_S421D_kif1a = polyfit(thresholded_green_all{2}, thresholded_red_all{2}, 1); %create linear regression
mdly_all_S421D_kif1a= polyval(mdl_all_S421D_kif1a,thresholded_green_all{2}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_S421D_kif1a = corr2(thresholded_green_all{2}, thresholded_red_all{2});
r_sqr_S421D_kif1a = power(correlation_coeff_S421D_kif1a,2);

mdl_all_WT_kif3a = polyfit(thresholded_green_all{3}, thresholded_red_all{3}, 1); %create linear regression
mdly_all_WT_kif3a= polyval(mdl_all_WT_kif3a,thresholded_green_all{3}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_WT_kif3a = corr2(thresholded_green_all{3}, thresholded_red_all{3});
r_sqr_WT_kif3a = power(correlation_coeff_WT_kif3a,2);

mdl_all_S421D_kif3a = polyfit(thresholded_green_all{4}, thresholded_red_all{4}, 1); %create linear regression
mdly_all_S421D_kif3a= polyval(mdl_all_S421D_kif3a,thresholded_green_all{4}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_S421D_kif3a = corr2(thresholded_green_all{4}, thresholded_red_all{4});
r_sqr_S421D_kif3a = power(correlation_coeff_S421D_kif3a,2);


mdl_all_WT_kif16b = polyfit(thresholded_green_all{5}, thresholded_red_all{5}, 1); %create linear regression
mdly_all_WT_kif16b= polyval(mdl_all_WT_kif16b,thresholded_green_all{5}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_WT_kif16b = corr2(thresholded_green_all{5}, thresholded_red_all{5});
r_sqr_WT_kif16b = power(correlation_coeff_WT_kif16b,2);

mdl_all_S421D_kif16b = polyfit(thresholded_green_all{6}, thresholded_red_all{6}, 1); %create linear regression
mdly_all_S421D_kif16b= polyval(mdl_all_S421D_kif16b,thresholded_green_all{6}); 
%this calculates the y values using the polynomial fit model
correlation_coeff_S421D_kif16b = corr2(thresholded_green_all{6}, thresholded_red_all{6});
r_sqr_S421D_kif16b = power(correlation_coeff_S421D_kif16b,2);

 figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{1}, thresholded_red_all{1});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{1}, mdly_all_WT_kif1a, 'color', 'b', 'LineWidth', 1); % plot the linear regression
                
figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{2}, thresholded_red_all{2});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{2}, mdly_all_S421D_kif1a, 'color', 'b', 'LineWidth', 1); % plot the linear regression
                
figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{3}, thresholded_red_all{3});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{3}, mdly_all_WT_kif3a, 'color', 'b', 'LineWidth', 1); % plot the linear regression
                
figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{4}, thresholded_red_all{4});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{4}, mdly_all_S421D_kif3a, 'color', 'b', 'LineWidth', 1); % plot the linear regression

figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{5}, thresholded_red_all{5});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{5}, mdly_all_WT_kif16b, 'color', 'b', 'LineWidth', 1); % plot the linear regression
                
figure
colormap(parula);
Nbins=100;
hold on, hb = binscatter(thresholded_green_all{6}, thresholded_red_all{6});
hb.NumBins = [Nbins Nbins];
hb.ShowEmptyBins = 'off';
hb.FaceAlpha = 0.75;
xlabel(xaxis); ylabel(yaxis);
set(gca,'FontSize',30);
pbaspect([1 1 1]);
colorbar;
hold on
% plot(thresholded_green_all{6}, mdly_all_S421D_kif16b, 'color', 'b', 'LineWidth', 1); % plot the linear regression
                

fileprefix='Motor_data_removedbad';

% T_WT_khc=table(model{1},r_squared{1});
% T_S421D_khc=table(model{2},r_squared{2});
% T_WT_dyn=table(model{3},r_squared{3});
% T_S421D_dyn=table(model{4},r_squared{4});

T_WT_kif1a=table(model{1},r_squared{1});
T_S421D_kif1a=table(model{2},r_squared{2});
T_WT_kif3a=table(model{3},r_squared{3});
T_S421D_kif3a=table(model{4},r_squared{4});
T_WT_kif16b=table(model{3},r_squared{3});
T_S421D_kif16b=table(model{4},r_squared{4});

% cd(save_dir);
% writetable(T_WT_khc,'WT_khc_fitting');
% writetable(T_S421D_khc,'S421D_khc_fitting');
% writetable(T_WT_dyn,'WT_dyn_fitting');
% writetable(T_S421D_dyn,'S421D_dyn_fitting');


% figure('Name','Number of Microtubule Associated Motors Violin','NumberTitle','off'), hold on,
% violin(r_squared, 'xlabel',{'WT khc','S421D khc','WT dyn','S421D dyn'},'facecolor',[0.5 0.5 0.5; 0 0.75 0.75; 0 0 0 ; 0 0.4 0.4],'edgecolor','k','mc','k','medc','k--');
% h1=notBoxPlot_nodotorline(r_squared{1},1,'style','line');
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(r_squared{2},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.75 0.75],'MarkerSize', 5);
% h3=notBoxPlot_nodotorline(r_squared{3},3,'style','line');
% set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0 0],'MarkerSize', 5);
% h4=notBoxPlot_nodotorline(r_squared{4},4,'style','line');
% set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
% publication_fig(0,0,1);
% ylim([-0.2 1.25]);
% pbaspect([2 1 1])
% ylabel('Correlation coefficient');
% 
% [WT_S421D_kin p_WT_S421D_kin]=ttest2(r_squared{1},r_squared{2})
% [WT_S421D_dyn p_WT_S421D_dyn]=ttest2(r_squared{3},r_squared{4})

figure('Name','Number of Microtubule Associated Motors Violin','NumberTitle','off'), hold on,
violin(r_squared, 'xlabel',{'WT kif1a','S421D kif1a','WT kif3a','S421D kif3a', 'WT kif16b', 'S421D kif16b'},'facecolor',[0.5 0.5 0.5; 0 0.75 0.75; 0 0 0 ; 0 0.4 0.4; 0.3 0.3 0.3; 0 0.55 0.55],'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(r_squared{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(r_squared{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.75 0.75],'MarkerSize', 5);
h3=notBoxPlot_nodotorline(r_squared{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0 0],'MarkerSize', 5);
h4=notBoxPlot_nodotorline(r_squared{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
h5=notBoxPlot_nodotorline(r_squared{5},5,'style','line');
set(h5.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize', 5);
h6=notBoxPlot_nodotorline(r_squared{6},6,'style','line');
set(h6.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.55 0.55],'MarkerSize', 5);
hLeg=legend('h1', 'h2', 'h3', 'h4', 'h5', 'h6');
set(hLeg,'visible','off');
publication_fig(0,0,1);
ylim([-0.5 1.25]);
pbaspect([2 1 1])
ylabel('Correlation coefficient');

[WT_S421D_kif1a p_WT_S421D_kif1a]=ttest2(r_squared{1},r_squared{2})
[WT_S421D_kif3a p_WT_S421D_kif3a]=ttest2(r_squared{3},r_squared{4})
[WT_S421D_kif16b p_WT_S421D_kif16b]=ttest2(r_squared{5},r_squared{6})

% figure('Name','Number of Microtubule Associated Dynein Violin','NumberTitle','off'), hold on,
% violin(r_squared{2,4}, 'xlabel',{'WT','S421D'},'facecolor',[0.5 0.5 0.5; 0 0.75 0.75],'edgecolor','k','bw',0.1,'mc','k','medc','k--');
% publication_fig(0,0,1);
% ylabel('Number of Microtubule-Associated Dynein');