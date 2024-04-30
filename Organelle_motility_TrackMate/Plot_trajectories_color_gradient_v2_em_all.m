%% Analyze Trajectories by color gradient

%Emily modifying to suit Muriel's data on Sept 17/22

% close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');

min_lr=0.6;
 
n_int=1;

% Create a color gradient matrix:
min_disp=0:n_int:13;
cmapp1=parula;
cmapp2=flipud(cmapp1); %Colormap goes from dark to light
gradient_matsz=size(cmapp2);
resiz_mat_1=round(gradient_matsz(1)./numel(min_disp));
resize_mat=1:resiz_mat_1:gradient_matsz(1);
cmapp5=[cmapp2(resize_mat,1),cmapp2(resize_mat,2),cmapp2(resize_mat,3)];
cmapp6=flipud(cmapp5);
cmapp6=[cmapp6];

% %Rab5 colours
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0 0.4 0.4];
% colour_WThtt=[0.3 0.3 0.3];
% colour_S421Dhtt=[0 0.4 0.68];

%Lyso colours
colour_WT=[0.729 0.831 0.957];
% colour_S421D=[0.757 0.867 0.776];
% colour_WThtt=[0.153 0.227 0.373];
% colour_S421Dhtt=[0.071 0.212 0.141];

tic
for k_choose = 1

% if k_choose == 1    % WT
    col1=colour_WT;
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_pos/');
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_pos/');
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_pos/');
%     cd('/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_pos/');
    fl='position.mat';
    an=0;
    DT=0.5;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    all_pos=[];
    
% end
save_fig_final='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Figures_motility/';

fls=dir(fl);
    
for k=1:numel(fls)
    display(fls(k).name)
    load(fls(k).name);
    %ab = ab(~cellfun(@isempty,{ab.position}));
    position=position(~cellfun(@isempty,(position)));
    for j=1:numel(position)
        kf=kf+1;
        range_disp(kf)=range(position{j});
        dat_pos_tim{kf}=position{j};
        %all_pos=[all_pos;ab(j).position];
    end
    
    for id=1:(numel(min_disp)+1)
        if id<min_disp(end)
        jk = find(range_disp > (min_disp(id)-(n_int./2)) & range_disp < (min_disp(id)+(n_int./2)));
        display(min_disp(id))
        display(min_disp(id)+(n_int./2))
        jk_ind{id}=jk;
        range_find{id}=range_disp(jk);
        data_position{id}=dat_pos_tim(jk);
        elseif id>min_disp(end)
        jk=find(range_disp > (min_disp(end)));
        jk_ind{id}=jk;
        range_find{id}=range_disp(jk);
        data_position{id}=dat_pos_tim(jk);
        end
    end
     
end
    
    dat_pos_tim(~cellfun('isempty',dat_pos_tim));
    data_position(cellfun('isempty', data_position)) = {{[0,0]}};
    all_pos=cell2mat(dat_pos_tim');

    out_pos_S421D_lyso=all_pos(find(all_pos>0));
    in_pos_S421D_lyso=all_pos(find(all_pos<0));
    frac_out_pos_S421D_lyso=numel(out_pos_S421D_lyso)/numel(all_pos);
    frac_in_pos_S421D_lyso=numel(in_pos_S421D_lyso)/numel(all_pos);
    
    %% Bin all the position points and plot them:
    pos_x=-20:0.05:20;
    pos_bin=hist(all_pos,pos_x);
    
    figure(2), hold on,
    bar(pos_x,pos_bin,'FaceColor',col1,'Facealpha',1,'BarWidth',0.5);
    set(gca, 'YScale', 'log');
    xlabel('Position (\mum)'); ylabel('Frequency'); 
    publication_fig(0,0,1);
    xlim([-15 15]);
    
    display(numel(dat_pos_tim))
    
    data_position=flipud(data_position');
%     
    for k2=1:numel(data_position)
%         hold on,
%         figure(1);
%         cellfun(@(x) plot(x, 'Color',cmapp6(k2,:)), data_position{k2});
%         publication_fig(0,0,1);
%         xlabel('Time (Frame)'); ylabel('Position (\mum)');
%         ylim([-20 20]); xlim([0 1500]);
%         hold on,
%         
%         hold on,
%         figure(3); %hold on,
%         subplot(1,2,1)
%         cellfun(@(x) plot(x, 'Color',cmapp6(k2,:)), data_position{k2});
%         publication_fig(0,0,1);
%         xlabel('Time (Frame)'); ylabel('Position (\mum)');
%         ylim([-20 20]); xlim([0 1500]);
%         hold on,
%         subplot(1,2,2)
%         bar(pos_x,pos_bin,'FaceColor',col1,'Facealpha',0.5,'BarWidth',0.5);
%         set(gca, 'YScale', 'log');
%         ylabel('Frequency'); 
%         publication_fig(0,0,1);
%         xlim([-20 20]);
%         view([90 -90]);

%         hold on,
%         figure(4); %hold on,
%         subplot(1,2,1)
%         cellfun(@(x) plot(x, 'Color',cmapp6(k2,:)), data_position{k2});
%         publication_fig(0,0,1);
%         xlabel('Time (Frame)'); ylabel('Position (\mum)');
%         ylim([-20 20]); xlim([0 1500]);
%         hold on,
%         subplot(1,2,2)
%         histogram(all_pos,'BinEdges',pos_x,'FaceColor',col1);
%         set(gca, 'YScale', 'log');
%         ylabel('Frequency'); 
%         publication_fig(0,0,1);
%         xlim([-20 20]);
%         view([90 -90]);
        
    end
    
%     fl_nm=['20221206_WT_r5_', 'all_Trajectories_color_coded'];
%     saveas(gca, fullfile(save_fig_final, fl_nm), 'svg');
%     saveas(gca, fullfile(save_fig_final, fl_nm), 'png');
%     saveas(gca, fullfile(save_fig_final, fl_nm), 'fig');
%     saveas(gca,fl_nm,'png')

    clear res all_pos pos_x pos_bin;
    clear pos tim dat_pos_tim xk_yk lr_p_2 lr_m_2 t_p_2 t_m_2 data_position range_find jk_ind dat_pos_tim;
    clear range_disp
    
    proc_runs=[];
    proc_times=[];
    nvproc=[];
    rng=[];
    all_pos=[];
%     clear fl_nm
    %close all
end

toc
