% upenn             23.01.2009                  adam g. hendricks
% plot and analyze force trace data
% Directions: 1.  run with fc=0 to write data for CPCpackage

%notes:
%100930 - calculate Beta based on frequency response at Fexc
%120120 - use calibration from trap_frf_r6.m
%       - calculate forces ant. / ret. from points off of image (mtoc,
%       bead)
clear all
close all
% addpath('~/Documents/postdoc/analysis/functions')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT/agh-optical-trap-codes-master 2/')
% disp('To Do: (1) apply std. dev. and other criteria - see Matt Lang paper
% (2) peak finder')

% options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zcont=1; %0-calibration, 1-stall force, 2-force feedback


% datdir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210731_OT_S421D_FN_late/'; %fs=20e3;
% datdir='/Volumes/Emily_htt_2/HEK_OT/20210926_WT_OT_early/'; %fs=20e3;
datdir='/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/'; %fs=20e3;
calib_fname={'.txt'};
image_fname={''};
force_fname={'20240416_kolf4_OT_s1_motors_5.txt'}; %
Beta= 62962.2232936428 ; %[nm/V] 
K_trap= 0.017908396407947  ; %[pN/nm] 
r_cam_ant=[-1 1]; %[dx dy] MT polarity on camera (pointing toward + end)
%
Fs_force=2e4; %Emily sampling frequency, actually not 2000 but 20000.
cht=0;
chxpd=2; % Emily make sure this references the correct column (reference Abdullah's code)
chypd=1; % Emily make sure this references the correct column (reference Abdullah's code)
chxaod=8;
% zcont=0; %Emily commented because I want zcont=1 to analyze force traces
% Gain=[10,10,10,10,1,1,2,1]; %gain on amplifier-Emily-Adam said to set these to 1
Gain=[1,1,1,1,1,1,1,1]; %gain on amplifier-Emily-Adam said to set these to 1
% flab='_traces_analyzed';%'_pLBC_force'; %append to data file

disp(['Directory: ',datdir])
disp(['Calibration filename: ',calib_fname])
disp(['Force trace filename: ',force_fname])

Nmdn=41; %window size for median filter
% Cstall_time=0.01;%0.1%0.01 %[s] minimum stall plateau time
Cstall_time=0.01;%0.1%0.01 %[s] minimum stall plateau time
Cstall_vel=50; %[nm/s] minimum stall velocity (up)
Cstall_snapvel=-100;%-500; %[nm/s] minimum stall velocity (dwn)

k_var=7; %variance scaling for plot

% constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp(['Beta = ', num2str(Beta), ' nm/V']);
    disp(['k_trap = ',num2str(K_trap),' pN/nm']);
    disp(['r_cam_ant = ', num2str(r_cam_ant), ' direction']);%Emily changed so you can actually see it

%load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdat=textread([datdir,force_fname{1}],'','headerlines',8);
[Nfdat,Mfdat]=size(fdat); %Nfdat is number of rows, Mfdat is number of columns.
ftime=[0:(Nfdat-1)]'./Fs_force; %[s] Emily-why do we divide the time by 20000? Oh the sampling frequency!

ftime_seg=buffer(ftime,120000);


% for ftime_seg_curr=1:10
%    
%     ftime=ftime_seg(:,ftime_seg_curr);
%Emily commented
% if cht==0
%     ftime = [0:(Nfdat-1)]'.*(1/Fs_force); %Emily-Isn't this the same as the previously set ftime?
% 
% else
%     ftime=Nfdat(:,cht);%
% end
%jkp=find(ftime<150); disp('keep only t<150 s')
%ftime=ftime(jkp);
%fdat=fdat(jkp,:);


xpd=Gain(chxpd).*fdat(:,chxpd);
ypd=Gain(chypd).*fdat(:,chypd);

xpd_seg=buffer(xpd,120000);
ypd_seg=buffer(ypd,120000);

%[xpdm,xpdvar]=medvarfilt1(xpd,Nmdn);
%[ypdm,ypdvar]=medvarfilt1(ypd,Nmdn);

file_orient=[force_fname{1}(1:end-4),'_orient.mat'];
% if exist(['./trap-calibration-data/',file_orient],'file')==2
%     z_orient=input('Use saved orientation data (0/1)? ');
%     if z_orient==1
%         load(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',file_orient]);
%     end
% else
%     z_orient=0;
% end
% if z_orient==0 %Emily- pretty sure this is the case for my data

    r_pd_ant=[r_cam_ant(1),r_cam_ant(2)].*[-1 -1]; %PD coordinate frame 
    r_pd_ant=r_pd_ant./sqrt(sum(r_pd_ant.^2)); %normalize
    r_pd_off=cross([r_pd_ant,0],[0 0 -1]);
    r_pd_off=r_pd_off(1:2); %off-axis vector

    %plot xpd, ypd signals - choose zero
    %load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fdat=textread([datdir,force_fname{1}],'','headerlines',8);
    [Nfdat,Mfdat]=size(fdat);
%     ftime=[0:(Nfdat-1)]'./Fs_force; %[s]

%     if cht==0 %Emily-is this not the same as on lines 102-109?
%         ftime=[0:(Nfdat-1)]'.*(1/Fs_force);
% 
%     else
%         ftime=fdat(:,cht);%
%     end
    %jkp=find(ftime<150); disp('keep only t<150 s')
    %ftime=ftime(jkp);
    %fdat=fdat(jkp,:);

    xpd=Gain(chxpd).*fdat(:,chxpd);
    ypd=Gain(chypd).*fdat(:,chypd);
    %[xpdm,xpdvar]=medvarfilt1(xpd,Nmdn);
    %[ypdm,ypdvar]=medvarfilt1(ypd,Nmdn);
    
%     for ftime_seg_curr=1:10
   
%     ftime=ftime_seg(:,ftime_seg_curr);
%     
%     xpd=xpd_seg(:,ftime_seg_curr);
%     ypd=ypd_seg(:,ftime_seg_curr);

    figure, subplot('position',[0.05 0.05 0.6 0.9]), plot(ftime,xpd)
    %hold on, plot(ftime,xpdvar,'color',[0.7 0.7 0.7])
    xlabel('time (s)'), ylabel('Xpd (V)')
    subplot('position',[0.75 0.05 0.2 0.9]),
    [n,x]=hist(xpd,100);
    barh(x,n)

    disp('Choose zero force region from histogram (Xpd) ')
    %[trefx,xref]=ginput(2);
    % jx0=find(ftime>trefx(1) & ftime<trefx(2));
    % xpd0=mean(xpd(jx0));
    [nref,xpd0]=ginput(1);


    figure, subplot('position',[0.05 0.05 0.6 0.9]), plot(ftime,ypd)
    %hold on, plot(ftime,ypdvar,'color',[0.7 0.7 0.7])
    xlabel('time (s)'), ylabel('Ypd (V)')
    subplot('position',[0.75 0.05 0.2 0.9]),
    [n,x]=hist(ypd,100);
    barh(x,n)
    disp('Choose zero force region from histogram (Ypd) ')
    %[trefy,yref]=ginput(2);
    % jy0=find(ftime>trefy(1) & ftime<trefy(2));
    % ypd0=mean(ypd(jy0));
    [nyref,ypd0]=ginput(1);

    rpd=[xpd-xpd0,ypd-ypd0]; %Emily-Normalizing to the user inputted zero force region
    
    %Adam's addition to colour the data by timepoint
    
%     Nseg = floor(Nfdat/10);
% 
%     figure, hold on
%     xlabel('Xpd'), ylabel('Ypd');
%     cmap = colormap(parula(10));
%     for kseg=1:9
%         jseg = kseg*Nseg + [1:Nseg];
%         plot(Beta.*rpd(jseg,1),Beta.*rpd(jseg,2),'.','color',cmap(kseg,:)), axis equal
%         %alpha(0.2)
%     end
    
    figure, plot(Beta.*rpd(:,1),Beta.*rpd(:,2),'.'), axis equal
    xlabel('Xpd'), ylabel('Ypd');
    disp('Choose axis of motility (start and end): ')
    [xmot,ymot]=ginput(2);
    rpt=[diff(xmot),diff(ymot)];
    rpt=rpt./sqrt(rpt(1).^2+rpt(2).^2);

    %ppt=polyfit(rpd(:,1),rpd(:,2),1);
    %rpt=[1 ppt(1)];
    if r_pd_ant*rpt' < 0 
        rpt=-rpt;
    end
    rpt_off=cross([rpt,0],[0 0 -1]);
    rpt_off=rpt_off(1:2);

    save(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',file_orient],'rpt','rpt_off','xpd0','ypd0')
% else
%     rpd=[xpd-xpd0,ypd-ypd0];
% end




%rpdc=[rpd*r_pd_ant',rpd*r_pd_off']; %on-axis (+ = ant.), off-axis (based on image)
rpdc=[rpd*rpt',rpd*rpt_off']; %on-axis (+ = ant.), off-axis (based on trajectory)

%Emily modifying to filter out data that leaves the trap (nonlinear regime)
%and to make the variables for the figures simpler
on_axis_pos= Beta.*rpdc(:,1);
off_axis_pos=Beta.*rpdc(:,2);
on_axis_force=Beta*K_trap.*rpdc(:,1);
off_axis_force=Beta*K_trap.*rpdc(:,2);
% for n=1:numel(on_axis_pos)
% if on_axis_pos(n)> 300
%     on_axis_pos(n)=NaN;
% else
% end

% end

% Emily- plotting on axis and off axis positions and forces 
figure(3), hold on, plot(on_axis_pos,off_axis_pos,'r.'), axis equal
xlabel('on-axis position (nm)'), ylabel('off-axis position (nm)')

figure(4), plot(ftime,on_axis_pos)
xlabel('time (s)'), ylabel('on-axis position (nm)')

figure(5), plot(ftime,off_axis_pos)
xlabel('time (s)'), ylabel('off-axis position (nm)')

figure(6), plot(ftime,on_axis_force);
xlabel('time (s)'), ylabel('on-axis force (pN)')

figure(7), plot(ftime,off_axis_force);
xlabel('time (s)'), ylabel('off-axis force (pN)')

if exist('kaod','var')==1 %Emily-This variable is not previously defined, so I think this condition does not apply
    xpdm=medfilt1(on_axis_pos,Nmdn,10e3);
    xaod=(kaod/Gain(chxaod)).*fdat(:,chxaod);
    xaodm=medfilt1(xaod,Nmdn,10e3);
    
    figure(100), plot(ftime,xaod,'markersize',6,'color',0.7.*[1 1 1])
    hold on, plot(ftime,xaodm,'-','color',[0.7 0 0],'linewidth',2)
    hold on, plot(ftime,+xaod,'markersize',6,'color',0.7.*[1 1 1])
    hold on, plot(ftime,xaodm+xpdm,'-','color',[0 0 0.7],'linewidth',2)
    xlabel('time (s)'), ylabel('laser position (red,nm) / bead position (blue,nm)')
    set(gca,'ytick',-1200:16:1200)
    set(gca,'ygrid','on')
end

zcont=input('Analyze force traces (0/1)? ');
if zcont==1 %analyze force trace dataf
    %[rpdcm,rpdcvar]=medvarfilt1(rpdc,Nmdn);     %median filter
    rpdcm=medfilt1(rpdc,Nmdn,10e3); %Emily-what is this factor of 10000?, block size
    %rpdcm=sgolayfilt(rpdc,3,Nmdn+1);

    %figure(3), hold on
    %plot(Beta.*rpdcm(:,1),Beta.*rpdcm(:,2),'c.')
    %scatter(Beta.*rpdcm(:,1),Beta.*rpdcm(:,2),[],ftime,'.')

    figure(4), hold on,
    plot(ftime,Beta.*rpdcm(:,1),'color','r','linewidth',1);

    figure(5), hold on,
    plot(ftime,Beta.*rpdcm(:,2),'color','r','linewidth',1);

    figure(6), hold on,
    plot(ftime,Beta*K_trap.*rpdcm(:,1),'color','r','linewidth',1);%[0.7 0.7 0.7])

    figure(7), hold on,
    plot(ftime,Beta*K_trap.*rpdcm(:,2),'color','r','linewidth',1);%[0.7 0.7 0.7])


    %     find stall events
    z_stall_events=1;%input('Find stall events (0/1)? ');
    if z_stall_events==1
        %disp('Set stall threshold.')
        %[xthreshold, Fthreshold]=ginput(1);
        Fthreshold=0.5; %Emily-the force threshold in pN.
        jstall=find(Beta*K_trap.*rpdcm(:,1)>Fthreshold | Beta*K_trap.*rpdcm(:,1)<-Fthreshold);
        jevent=find(diff(jstall)>10); %Emily-Not sure I get this exactly
        jevent=[0;jevent;length(jstall)];
        for kevent=2:length(jevent)
            ji=jstall(jevent(kevent-1)+1);
            jf=jstall(jevent(kevent));

            %         hold on, plot(ftm(ji),ffm(ji),'c*',ftm(jf),ffm(jf),'co','MarkerSize',20)
            stall_ti(kevent-1)=ftime(ji);
            stall_tf(kevent-1)=ftime(jf);
            if mean(Beta*K_trap.*rpdcm(ji:jf,1))>0
                stall_force(kevent-1)=max(Beta*K_trap.*rpdcm(ji:jf,1));
                %figure(6), hold on, plot(ftime(ji:jf),Fthreshold.*ones(size(ftime(ji:jf))),'c.')
            else
                stall_force(kevent-1)=min(Beta*K_trap.*rpdcm(ji:jf,1));
                %figure(6), hold on, plot(ftime(ji:jf),-Fthreshold.*ones(size(ftime(ji:jf))),'c.')
            end
            stall_time(kevent-1)=ftime(jf)-ftime(ji); %plateau time
            jvel=find(ftime>(ftime(ji)-Cstall_time) & ftime<(ftime(ji)+Cstall_time));% Emily- I think this is finding the time of the stall event
            pvel=polyfit(ftime(jvel),Beta.*rpdcm(jvel,1),1);
            stall_vel(kevent-1)=pvel(1); % stall velocity
            jsnap=find(ftime>(ftime(jf)-0.05) & ftime<(ftime(jf)+0.05));
            psnap=polyfit(ftime(jsnap),Beta.*rpdcm(jsnap,1),1);
            stall_snapvel(kevent-1)=psnap(1); %snap back velocity
        end

        %     apply conditions
        jkeep=find(stall_time > Cstall_time & abs(stall_vel) > Cstall_vel & -abs(stall_snapvel) < Cstall_snapvel);
%         for k=1:length(jkeep),
%             if mean(stall_force(jkeep(k))>0),
%                 figure(6), hold on, plot([stall_ti(jkeep(k)),stall_tf(jkeep(k))],[Fthreshold,Fthreshold],'m-','LineWidth',2)
%             else
%                 figure(6), hold on, plot([stall_ti(jkeep(k)),stall_tf(jkeep(k))],-[Fthreshold,Fthreshold],'m-','LineWidth',2)
%             end
%         end
        if numel(jkeep)>3
            figure, hist(stall_force(jkeep),round(numel(jkeep)/3)), xlabel('stall force (pN)'), ylabel('Count')
        end

        stall_force_keep=stall_force(jkeep);
%%

        % pairwise distances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         zpw=input('Calculate pairwise distances (0/1)? ');
%         if zpw==1,
%             Fthresh_pairwise=5;
%             jskip=10;
%             Nbins=100;
% 
%             %use same data that was previously fit using Kerssemakers algorithm:
%             %Emily-what is this data?
%             flist_steps=dir(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',force_fname{1}(6:end-4),'*_steps.mat']);
% 
%             % sum histograms, calc. all pairs
%             xbins=-200:0.2:200;
%             for ks=1:numel(flist_steps),
%                 disp(['pairwise calculation - ',num2str(ks),'/',num2str(numel(flist_steps))])
%                 dats=load(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',flist_steps(ks).name]);
%                 jpw=dats.jstep;
%                 %xpw=Beta.*sqrt(rpdcm(jpw,1).^2 + rpdcm(jpw,2).^2);
%                 xpw=Beta.*rpdcm(jpw,1);%interp1(ftime(jpw),rpdcm(jpw,1),ftime(jpw(1:jskip:end)),'pchip');
%                 Nx=numel(xpw);
%                 npw=zeros(size(xbins));
%                 for kxi=1:Nx,
%                     jk=[(kxi+1):Nx];
%                     dpwk=xpw(jk)-xpw(kxi);
% 
%                     npwk=hist(dpwk,xbins);
%                     npw=npw+npwk;
%                     clear dpwk
%                 end
% 
%                 figure, plot(xbins,npw,'.-')
%                 set(gca,'xlim',[-32 32])
%                 set(gca,'xtick',[-32:4:32]), grid on
%                 title('Pairwise Displacement')
% 
%                 save(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',flist_steps(ks).name(1:end-4),'_pw.mat'],'npw','-append');
% 
%                 clear dats jpw npw
%             end
% 
%             % sum distances, calc. every jskip-th pair
%             %         for ks=1:numel(flist_steps),
%             %             disp(['pairwise calculation - ',num2str(ks),'/',num2str(numel(flist_steps))])
%             %             dats=load(['./trap-calibration-data/',flist_steps(ks).name]);
%             %             jpw=dats.jstep;
%             %             %xpw=Beta.*sqrt(rpdcm(jpw,1).^2 + rpdcm(jpw,2).^2);
%             %             xpw=Beta.*interp1(ftime(jpw),rpdcm(jpw,1),ftime(jpw(1:jskip:end)),'pchip');
%             %             Nx=numel(xpw);
%             %             dpw=[];
%             %             for kxi=1:Nx,
%             %                 jk=[(kxi+1):Nx];
%             %                 dpwk=xpw(jk)-xpw(kxi);
%             %                 dpw=[dpw;dpwk];
%             %             end
%             %
%             %             xbins=[-32:0.2:32];
%             %             n_pw=hist(dpw,xbins);
%             %
%             %             figure, plot(xbins,n_pw,'.-')
%             %             set(gca,'xtick',[-32:4:32]), grid on
%             %             title('Pairwise Displacement')
%             %
%             %             save(['./trap-calibration-data/',flist_steps(ks).name(1:end-4),'_pw.mat'],'dpw');
%             %
%             %             clear dats jpw dpw n_pw
%             %         end
% 
%             %         figure, bar(xbins,n_pw_ret,1)
%             %         set(gca,'xtick',[-28:4:28]), grid on
%             %         title('Retrograde Pairwise Distance')
%         end
% %%
%         
%         
% %%
%         % gaussian kernel density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         zgk=input('Calculate gaussian kernel density (0/1)? '); 
%         if zgk==1,
%             gk_win=3; %Emily-Gaussian kernel window size 
%             Fx=1/0.2; %[1/nm] Emily, what is this constant?
% 
%             %use same data that was previously fit using Kerssemakers algorithm:
%             %Emily-do I have this data? I don't think so.
%             flist_steps=dir(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',force_fname{1}(6:end-4),'*_steps.mat'])
%   
%             for ks=1:numel(flist_steps),
%                 disp(['gaussian kernel calculation - ',num2str(ks),'/',num2str(numel(flist_steps))])
%                 dats=load(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',flist_steps(ks).name])
% 
%                 % gaussian kernel filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                 xdatk=Beta.*(rpdcm(dats.jstep,1)-rpdcm(dats.jstep(1),1));
%                 Ngk=numel(dats.jstep);
%                 Nk=floor(Ngk/gk_win);
%                 sx=sign(xdatk(end)); %Emily-gives the sign (+,-, 0) of the rpdcm data.
%                 xgk=0:sx*(1/Fx):xdatk(end);
%                 ngk=zeros(size(xgk));
%                 for kg=1:Nk,
%                     jk=(kg-1)*gk_win + [1:gk_win];
%                     [muk,sigk]=normfit(xdatk(jk));
%                     %tk=mean(jk)/Fs;
%                     ngkk=normpdf(xgk,muk,sigk); %Emily-what is the purpose of this whole gk function?
%                     if sum(isnan(ngkk))==0,
%                         ngk=ngk+ngkk;
%                     end
%                 end
%                 
%                 Nfft=2^nextpow2(numel(xgk));%80*5;%numel(jk);
%                 Nwindow=min([8*50,numel(xgk)]);
%                 Noverlap=floor(Nwindow/2);
%                 Nint=3;
%                 [ps_gk,psd_gk]=power_spectrum([xgk;ngk],Nwindow,Noverlap,Nfft,Nint,Fx);
%                 [c_gk,lags_gk]=xcorr(ngk);
% 
% %                 figure,
% %                 subplot(121), plot(ftime(dats.jstep),xdatk,'linewidth',2)
% %                 set(gca,'ylim',sx.*[0 100])
% %                 set(gca,'ytick',sx.*[0:8:100])
% %                 set(gca,'ygrid','on')
% %                 xlabel('Time (s)'), ylabel('Position (nm)')
% %                 subplot(122), plot(ngk./sum(ngk),xgk,'linewidth',2)
% %                 set(gca,'ylim',sx.*[0 100])
% %                 set(gca,'ytick',sx*[0:8:100])
% %                 set(gca,'ygrid','on')
% %                 xlabel('Position (nm)'), ylabel('Fraction')%, title('Gaussian Kernel Filter')
% % 
% %                 figure, loglog(ps_gk(:,1),ps_gk(:,2))
% %                 set(gca,'xtick',1./[16:-4:4,2,1])
% %                 set(gca,'xgrid','on')
% %                 set(gca,'xticklabel',[16:-4:4,2,1])
% %                 set(gca,'xlim',1./[32 1])
% %                 xlabel('$f_x (nm^{-1})$'), ylabel('Power Spectrum')
% 
% %                 figure, plot(lags_gk./Fx,c_gk)
% %                 xlabel('$\Delta x~$ (nm)'), ylabel('Autocorrelation')
% %                 set(gca,'xlim',[-32,32])
% %                 set(gca,'xtick',[-32:8:32])
% %                 set(gca,'xgrid','on')
% 
%                 save(['/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/HEK_OT/20210704_OT_WT_fibronectin/',flist_steps(ks).name(1:end-4),'_pw.mat'],'xgk','ngk','ps_gk','psd_gk','-append');
% 
%                 clear dats xgk ngk ps_gk psd_gk kg
%             end
%         end
% %%
% zfcyt=input('Calculate Fcyt (0/1)? '),
% if zfcyt==1,
% %         %force-velocity curves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         N=3; %order of polynomial fit
%         F=501; %window length
%         [SG0,SG1]=sgolayfiltN(N,F,1e4,1/Fs_force,rpdc(:,1));
% %         [SG0,SG1]=sgolayfilt(rpdc(:,1),N,F,1e4,1/Fs_force);%Emily changed, what is Fs_force? (currently 2e3)
% %         SG0=sgolayfilt(rpdc(:,1),N,F,1e4,1/Fs_force);%Emily changed, what is Fs_force? (currently 2e3)
%         
%         figure
%         subplot(2,1,1);
%         plot(ftime,, 'b',ftime,Beta.*SG0','r')
%         xlabel('time (s)'), ylabel('position')
% 
%         subplot(2, 1, 2);
%         plot(ftime,[Beta.*SG1'])
%         xlabel('time (s)'), ylabel('velocity')
%   
% %         F0=K_trap*Beta.*SG0;
% %         figure, plot(F0,Beta.*SG1,'.')
% %         xlabel('force'), ylabel('velocity')
% % 
% %         df=1;
% %         ff=-30:df:30;
% %          jpos=find(sign(SG0)==sign(SG1))
% %         for kf=1:numel(ff),
% %             jk=find(F0(jpos)>(ff(kf)-df/2) & F0(jpos)<(ff(kf)+df/2));
% %             vv(kf)=mean(SG1(jpos(jk)));
% %         end
% % 
% %         figure, plot(ff,vv,'.-')
% %         xlabel('force'), ylabel('velocity')
% %         
% %     % pure viscosity
% %     Fgamma=cdat.gammar*cdat.ig.Gamma0*Beta.*SG1;
% %     Fk0=cdat.k_cyt0*Beta.*SG0;
% %     Fk1=cdat.k_cyt1.*((Beta.*SG1).^(cdat.alpha)).*((Beta.*SG0).^(1-cdat.alpha)); % ~k*w^alpha
%     
% %     w0=2*pi.*100;
% %     G0=(1/6*pi*cdat.Rbead).*(cdat.k_cyt0 + cdat.k_cyt1.*((j.*w0).^cdat.alpha)./gamma(cdat.alpha)...
% %     + cdat.gammar*cdat.ig.Gamma0.*j.*w0);
% %     v0=imag(G0)./w0;
%     
% %     F_viscous=v0*Beta.*SG1;
%     %s=tf('s');
%     %Fcyt2dx=cdat.gammar*cdat.ig.Gamma0*s+cdat.k_cyt0+cdat.k_cyt1*s^cdat.alpha/gamma(cdat.alpha);
%     
%     %approximate transfer function of Fcyt/dx (integer powers of s)
% %%
%     ww=2*pi.*[0.1:0.1:5e3];%cdat.Fexc; Emily-what is this constant?
%     ff=ww./(2*pi);
%     
% %     Hgamma=tf([cdat.gammar*cdat.ig.Gamma0 0],[1]);
% %     Fgamma=lsim(Hgamma,Beta.*rpdcm(:,1),ftime,0);
% %     
% %     Hkc0=tf([cdat.k_cyt0],[1]);
% %     Fkc0=lsim(Hkc0,Beta.*rpdcm(:,1),ftime,0);
% %     
% %     hh_kcyt1=cdat.k_cyt1.*((j.*ww).^cdat.alpha)/(gamma(cdat.alpha));
% %     [b_kc1,a_kc1]=invfreqs(hh_kc1,ww,10,10);
% %     Hkc1=tf(b_kc1,a_kc1);
% %     Fkc1=lsim(Hkc1,Beta.*rpdcm(:,1),ftime,0);
%      
%     hh=cdat.gammar*cdat.ig.Gamma0*j.*ww+cdat.k_cyt0+cdat.k_cyt1.*((j.*ww).^cdat.alpha)/(gamma(cdat.alpha));
%     %hh=cdat.gammar*cdat.ig.Gamma0*j.*ww; %gamma only
%     %hh=cdat.k_cyt0+cdat.k_cyt1.*((j.*ww).^cdat.alpha)/(gamma(cdat.alpha)); %kcyt only 
%     [b,a]=invfreqs(hh,ww,8,8,[],10,0.01);%,'trace');
%     
%     Hcyt=tf(b,a);
%     tsim=0:1e-5:ftime(end);
%     Fcyt=lsim(Hcyt,interp1(ftime,Beta.*SG0,tsim),tsim,0);
%     Fcyt=interp1(tsim,Fcyt,ftime);
%     
% %     hh=(1./(j.*ww)).*hh; % in terms of velocity F/V
% %     [b,a]=invfreqs(hh,ww,4,4,[],1e3);
% %     Hcyt=tf(b,a);
% %     Fcyt=lsim(Hcyt,Beta.*SG1,ftime,0);
%     
%     %Fanalytic=cdat.gammar*cdat.ig.Gamma0*Beta.*SG1; %Fg=gamma*velocity
% %     Fanalytic=cdat.gammar*cdat.ig.Gamma0*Beta.*SG1...
% %         + cdat.k_cyt0*Beta.*SG0...
% %         + (cdat.k_cyt1/gamma(cdat.alpha)).*(SG1./SG0).^cdat.alpha.*(Beta.*SG0);
%     jan=find(ftime>=175 & ftime<=200);
%     Nan=numel(jan);
%     tan=ftime([1:Nan]);
%     H_kcyt0=zeros(size(tan));
%     H_kcyt0(1)=cdat.k_cyt0;
%     H_kcyt1=cdat.k_cyt1.*tan.^(-(cdat.alpha-1));
%     %H_kcyt1=H_kcyt1./H_kcyt1(2);
%     %H_kcyt1(1)=1;
%     F_kcyt0=conv(H_kcyt0,Beta.*SG0(jan));
%     F_kcyt1=conv(H_kcyt1,Beta.*SG0(jan));
%     F_gamma=cdat.gammar*cdat.ig.Gamma0*Beta.*SG1(jan);
%     Fanalytic=F_kcyt0(1:Nan)+F_kcyt1(1:Nan)+F_gamma;
%     
%     g=cdat.k_cyt0 + cdat.k_cyt1.*((j.*ww).^cdat.alpha)./gamma(cdat.alpha)...
%         + cdat.gammar*cdat.ig.Gamma0.*j.*ww;
%     Gp=g./(6*pi.*cdat.Rbead);
%     v=imag(g)./ww;
%     
%     figure, plot(ff,v,'k','linewidth',2)
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     
%     figure, gv1=gcf; 
%     subplot(311), hold on, ylabel('$F_{tot}$ (pN)')
%     subplot(312), hold on, ylabel('Velocity (nm/s)')
%     subplot(313), hold on, xlabel('Time (s)'), ylabel('$F_{\gamma}$ (pN)')
%     fv=[1,10,100];
%     cmapg=colormap(lines(numel(fv)));
%     Ft=K_trap*Beta.*SG0(jan);
%     figure(gv1), subplot(311), plot(ftime(jan),Ft,'k','linewidth',2)
%     figure(gv1), subplot(312), plot(ftime(jan),Beta.*SG1(jan),'k','linewidth',2)
%     for kv=1:numel(fv),
%        vk=interp1(ff,v,2*pi*fv(kv));
%        Fvk=vk*Beta.*SG1(jan);
%        
%        figure(gv1), subplot(311), hold on
%        plot(ftime(jan),Ft+Fvk,'color',cmapg(kv,:),'linewidth',1)
%        subplot(313), hold on
%        plot(ftime(jan),Fvk,'color',cmapg(kv,:),'linewidth',2)
%     end
%         
%     %Hcyt=ss(Hcyt); %put in state space form
%     %%Hcyt=frd(hh,ww);
%     %Hdat=ssdata(Hcyt);
%     %Ns=numel(Hdat(:,1));
%     %Fcyt=lsim(Hcyt,Beta.*rpdcm(:,1),ftime,zeros(Ns,1));
%  
%     [mtot,ptot]=bode(Hcyt,ww);    
%     figure, subplot(211), hold on,
%     semilogx(ff,20.*log10(abs(hh)),'b','linewidth',2)
%     semilogx(ff,20.*log10(mtot(:)),'c','linewidth',2)
%     subplot(212), hold on,
%     semilogx(ff,(180/pi).*angle(hh),'b','linewidth',2)
%     semilogx(ff,mod((180/pi).*ptot(:),360),'c','linewidth',2)
%     
% 
%      
% %     figure, 
% % %     subplot(411), plot(ftime,Fgamma), ylabel('$F_{\gamma}$')
% % %     subplot(412), plot(ftime,Fkc0), ylabel('$F_{kcyt,0}$')
% % %     subplot(413), plot(ftime,Fkc1), ylabel('$F_{kcyt,1}$')
% % %     subplot(414), 
%     figure, plot(ftime,Beta*K_trap.*SG0,'b',ftime, Beta*K_trap.*SG0+Fcyt','c'), xlabel('time (s)'), ylabel('$F_{trap}+F_{cyt}~$ (pN)')
%     hold on, plot(ftime(jan),Beta*K_trap.*SG0(jan)+Fanalytic,'g')
%     legend('$F_{trap}$','$F_{tot,invfreqs}$','$F_{tot,analytical}$')
%     
%         % Numerically calculate impulse response using ifft.m
%         jfc=jan;%find(ftime>=50 & ftime<=60);
%         Nx=numel(jfc);
%         ww=2*pi.*(0.5*Fs_force*linspace(0,1,Nx/2));
%         %hh=-(cdat.gammar*cdat.ig.Gamma0*j.*ww + cdat.k_cyt0);%+cdat.k_cyt1.*((j.*ww).^cdat.alpha)/(gamma(cdat.alpha)));
%         hh=-(cdat.gammar*cdat.ig.Gamma0*j.*ww); %viscosity only 
%         ht=Fs_force.*ifft(hh,Nx,'symmetric');
%         Fcyt2=(1/(Fs_force)).*conv(real(ht),Beta.*SG0(jfc));
%         Fcyt2=Fcyt2(1:Nx);
%         
%         hold on, plot(ftime(jfc),Beta*K_trap.*SG0(jfc)+Fcyt2,'m')%, xlabel('time (s)'), ylabel('$F_{cyt}~$ (pN)')
%          %legend('$F_{trap}$','$F_{tot,invfreqs}$','$F_{tot,\gamma \cdot v}$','$F_{tot,ifft}$')
% 
%          %in terms of velocity 
% %         hhv=(1./(j.*ww)).*hh;%-(cdat.gammar*cdat.ig.Gamma0.*ones(size(ww)));
% %         hvt=ifft(hh,Nx,'symmetric');
% %         Fcytv=conv(real(hvt),Beta.*SG1(jfc));
% %         Fcytv=Fcytv(1:Nx);
% %         
% %        hold on, plot(ftime(jfc),Fcytv,'r',ftime(jfc),cdat.gammar*cdat.ig.Gamma0*Beta.*SG1(jfc),'c')
% end
% %%
%         
%         
%         %plot combined force/position plot
%         figure(10), [ax,h1,h2]=plotyy(ftime,,ftime,Beta*K_trap.*rpdcm(:,1))
% 
%         set(ax(1),'ylim',(1/K_trap).*[-25 25])
%         %set(ax(1),'ytick',[-600:8:600])
%         %set(ax(1),'ygrid','on')
%         %set(ax(1),'ticklength',[0 0])
%         %set(ax(1),'ylabel','Position (nm)')
%         %set(ax(1),'xlabel','Time (s)')
% 
%         set(ax(2),'ylim',[-25 25])
%         %set(ax(2),'ytick',[-24:2:24])
%         %set(ax(2),'ylabel','Force (pN)')
%         color_force_traces_em(gcf); % Emily modified this function slightly
% 
%         z_step=input('Find steps (0/1)? ')
%         if z_step==1,
%             %call stepfinder
%             addpath('./smtracking')
%             figure(10)
%             disp('Choose time interval for step finding: ')
%             tstep=ginput(2);
%             jstep=find(ftime>tstep(1,1) & ftime<tstep(2,1));
%             jjstep=jstep;%(1:Nmdn:end);
%             [steps.fit,steps.steps,steps.histofsteps,steps.All_Steppedness]=stepfinder_force_trace_r6([ftime(jjstep),Beta.*rpdcm(jjstep,1)]);%Need stepfinder function Emily
% 
%             figure, plot(ftime(jstep),Beta.*rpdcm(jstep,1),'color',[0 0.5 0],'linewidth',2)
%             hold on, plot(ftime(jjstep),steps.fit,'color',[0.5 0 0],'linewidth',2)
%             set(gca,'ytick',[-200:8:200])
%             dx=[-24:2:24];
%             ndx=hist(steps.steps(5,:),dx);
%             figure, bar(dx,ndx,1), xlabel('Step Size (nm)'), ylabel('Counts')
%             set(gca,'xtick',[-24:4:24])
%         end

        zs=input('Save (0/1)? ');
        if zs==1
            date_now = datestr(now, 'dd.mm.yy-HH.MM.SS');
            fnsave=[date_now,force_fname{1}(1:end-4),'.mat'];
            save(['/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/',fnsave],'stall_force','jkeep','stall_time','stall_vel','stall_snapvel','stall_ti','stall_tf','ftime','on_axis_force','off_axis_force','K_trap','on_axis_pos','off_axis_pos')
        end
    end

    %     figure
    %     [ax,h1,h2]=plotyy(ftime,fpos,ftime,ff);
    %     ymax=max(fpos)+20;
    %     ymin=min(fpos)-20;
    %     set(ax(1),'YLim',[ymin ymax]);
    %     set(get(ax(1),'YLabel'),'String','Position (nm)')
    %     set(ax(1),'Ytick',[round(ymin):8:ymax]), grid on
    %     set(ax(2),'YLim',k_trap.*[ymin ymax])
    %     set(ax(2),'Ytick',[-100:10:100])
    %     set(get(ax(2),'YLabel'),'String','Force (pN)')



    %     %     to plot all stall force data
    %     dd=dir([datdir,'*stall_force.txt']);
    %     %     jp=7:14;
    %     FS=[];
    %     for k=1:length(dd),
    %         dk=textread([datdir,dd(k).name]);
    %         FS=[FS,dk(1:(end-1))];
    %     end

    %     [muhat,sigmahat]=normfit(FS);
    %     ffit=(min(FS)-0.5):0.1:(max(FS)+0.5);
    %     nfit=normpdf(ffit,muhat,sigmahat);
    %     [Ndat,Fdat]=hist(FS);
    %     figure, bar(Fdat,Ndat./sum(Ndat))
    %     hold on, plot(ffit,0.7.*nfit,'r-','LineWidth',2)
    %     xlabel('Stall Force (pN)'), ylabel('probability density')
    %     %     disp('Choose location for text on plot.')
    %     text(mean(ffit)+1.5, 0.3,{['N = ',num2str(sum(Ndat))],['\mu = ', num2str(muhat)],['\sigma = ',num2str(sigmahat)]})
end
%     end
% end