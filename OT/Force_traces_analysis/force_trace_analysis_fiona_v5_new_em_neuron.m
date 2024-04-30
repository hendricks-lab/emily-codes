% upenn             23.01.2009                  adam g. hendricks
% plot and analyze force trace data
% Directions: 1.  run with fc=0 to write data for CPCpackage

%notes:
%100930 - calculate Beta based on frequency response at Fexc

clear all, close all, clc
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/OT');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Abdullah_codes/Abdullah-codes-master/Optical Trapping/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
set(0,'DefaultFigureWindowStyle','docked')
% disp('To Do: (1) apply std. dev. and other criteria - see Matt Lang paper (2) peak finder')

% options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zcont=1; %0-calibration, 1-stall force, 2-force feedback

% addpath('/Volumes/Emily_htt_2/Neuron/isoHD18Q/20230209_18Q_lyso_bdnf_OT/OT/');
% datdir='/Volumes/Emily_htt_2/Neuron/isoHD18Q/20230209_18Q_lyso_bdnf_OT/OT/';
addpath('/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/');
datdir='/Volumes/Emily_2022/Neurons/20240416_kolf4_OT/';
% addpath('/Volumes/Emily_2022/Omar_OT/U2OS/20240418_U2OS_OT/');
% datdir='/Volumes/Emily_2022/Omar_OT/U2OS/20240418_U2OS_OT/';
% addpath('/Volumes/Emily_htt_2/HEK_OT/20210813_WT_OT_early_late');
% datdir='/Volumes/Emily_htt_2/HEK_OT/20210813_WT_OT_early_late/'; 
% addpath=('/Volumes/Emily_htt_2/Neuron/isoHD45Q/20220213_30Q_DIV7_OT/');
% datdir='/Volumes/Emily_htt_2/Neuron/isoHD30Q/20220213_30Q_DIV7_OT/';
%cd('C:\Users\abdul\Desktop\All the data - LBC\optical trap data\160928-abdullah\test\');
calib_fname={'none.txt'};
force_fname={'20240416_kolf4_OT_s1_motors_5.txt'};
%MT=xlsread('MT1.xlsx');

%mt_coords=[MT(:,6),MT(:,7)];
%mtx=mt_coords(:,1);
%mty=mt_coords(:,2);
%mty=-1*mt_coords(:,2);
zdim='xy';

% if range(mtx)>range(mty)
%    zdim='x';
% elseif range(mty)>range(mtx)
%     zdim='y';
% else
%     zdim='xy';
% end

Fs=20e3; %calibration data
Fs_force=20e3;   % Force measurements, Emily force measurement accidentally at 20e3, modified from 2e3 (what it should have been)
%PD_sensitivity=2.472e3; %calib01-1: wt01-1:7
%k_trap=0.00786;
%PD_sensitivity = 2412.2; %cal01-3: wt01-8:12
%k_trap = 0.00723;
%PD_sensitivity = 2473.0; %mean(cal02-1): sma01-1:3
%k_trap = 0.00786;
PD_sensitivity = 62962.2232936428; % Emily calculated using Loic's code %mean(cal02-2,4): sma01-6:14, for neuron data this is wrong but don't have good calibration
k_trap = 0.017908396407947; %Emily calculated using Loic's code, for neuron data this is wrong but don't have good calibration
k_trap1=k_trap;

cht=0;
chxpd=1;
Gxpd=1;
chypd=2;
Gypd=1;
chxaod=5;
Gaod=1;
chxstage=4;
Gxstage=1;

Fexc=0; %[Hz] frequency of excitation
%Fs=50e3;%20e3; %for calibration data
Rbead=0.1e3; %[nm] bead radius
%falias=15e3;
falias=4.8e3; %[Hz] anti-aliasing filter (correction to Lorentzian fit)
Qo=1e-5;%
f0=0;%1;%20;%5; %[Hz] building vibrations (0 if none)
fmin=5e1; % range to start fit
fmax=5e3;
Fexc_window=max(20,floor(Fexc/10));%150;%%250;%15;

k_exc=20.37*(16/16); %[nm/V] AOD, for frequency multiplier = 16
%chxst=5; disp('Laser Excitation')
% k_exc=10e3; %[nm/V] excitation scale factor
% chxst=4; disp('Stage Excitation')

disp(['Directory: ',datdir])
disp(['Calibration filename: ',calib_fname])
disp(['Force trace filename: ',force_fname])

% datdir='../data/fiona-data/100310-kinesin/';
% calib_fname='fiona100310_037.txt';
% force_fname='fiona100310_035.txt';
%k_trap=0.609;%1.07;%0.0374;%1.12;%0.0374; %[pN/nm] - enter 0 if you need to run calibration
%PD_sensitivity=73.3;%28.7;%284;%68.3;%284;%750; %[nm/V] photodiode sensitivity (initial guess)
% 20032009 - 525, 11032009 - 959.7, ??? - 419.4, 21052009 - 468.3D
%Dfit=10.11; %[V^2/s] (from start_fit.m)
% zrun=input('Input: calibration (0) or force (1)? ');
Nmdn=30; %window size for median filter
Cstall_time=0.01;%0.1%0.01 %[s] minimum stall plateau time
Cstall_vel=10; %[nm/s] minimum stall velocity (up)
Cstall_snapvel=10;%-500; %[nm/s] minimum stall velocity (dwn)

k_var=7; %variance scaling for plot

Nwindow=512*floor(Fs/2e3);%512*12;%512*6;floor(Fs/10);%
Noverlap=ceil(Nwindow/2);
Nfft=Fs;
Nint=10;

% Nwindow=512*12;%512*6;
% Noverlap=ceil(Nwindow/2);
% Nfft=20e3;
% Nint=10;
% Fs=20e3; %for calibration data


% constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kT=4.1; %1.381e-23*(25+273.15)*(1e12*1e9); %[pN*nm] (assuming 25 degC)
viscosity=1e-9;%1e-3*(1e12)/(1e9)^2; %[pN*s/nm^2] viscosity of water

% analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if zcont==0;%k_trap==0,  % put calibration data in file to be read by CPCpackage/start_fit.m
    %caldat=textread([datdir,calib_fname]);

    tpsx=[];
    psx=[];
    tpsdx=[];
    psdx=[];

    tpss=[];
    pss=[];
    tpsds=[];
    psds=[];
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Data/');
    for kd=1:numel(calib_fname);
        caldat=textread(calib_fname{kd},'','headerlines',8);
        Ncdat=length(caldat(:,chxpd));
        Tdrive=1/Fexc;
        %NTdrive=floor(caldat(end,chxpd)/Tdrive);
        
        ik=1:numel(caldat(:,chxpd));%80e3;%find(caldat(:,1)<=Tmsr);
        if cht==0,
            ctime=ik'./Fs;
        else
        ctime=caldat(ik,cht);
        end
        NTdrive=ctime(end)/Tdrive;
        Tmsr=Tdrive*NTdrive;
        cposx=Gxpd.*caldat(ik,chxpd); 
        cposy=Gypd.*caldat(ik,chypd);
        if Fexc~=0,
        cstagex=Gxstage.*caldat(ik,chxst);
        end
        %Nfft=floor(numel(ctime));

        %Noverlap=0;%floor(Nfft/2);
        %         A2=.01;
        %         f2=15;
        %         cexc=A2*sin(2*pi*f2*ctime); disp('DEBUGGING - see line 51')
        %         cposx=cposx+cexc;

        %     % power spectrum
        %     Nfft=[];%2^12;
        %     Nwindow=[];%1e3;
        %     Noverlap=[];%floor(Nwindow/2);
        %     Nint=10;
        %     [cps,cpsd]=power_spectrum([ctime,cpos]',Nwindow,Noverlap,Nfft,Nint,cfs);
        %     figure, loglog(cpsd(:,1),cpsd(:,2),'b.')
        %
        %     disp('Click on corner frequency')
        %     [fc,yc]=ginput(1);
        %     hold on, loglog(fc,yc,'ro')
        %     k_trap=2*pi*gamma_bead*fc;
        %     disp(['k_trap = ',num2str(k_trap),' pN/nm'])
        disp('Calculating power spectra ...')

        [psx_k,psdx_k]=power_spectrum([ctime,cposx]',Nwindow,Noverlap,Nfft,Nint,Fs);
        Tmsr=ctime(end);%caldat(end,1);
        
        tpsx=[tpsx,psx_k(:,1)];
        psx=[psx,psx_k(:,2)];
        tpsdx=[tpsdx,psdx_k(:,1)];
        psdx=[psdx,psdx_k(:,2)];
        Adrive1=0;
        if Fexc~=0
            [pss_k,psds_k]=power_spectrum([ctime,cstagex]',Nwindow,Noverlap,Nfft,Nint,Fs);

            tpss=[tpss,pss_k(:,1)];
            pss=[pss,pss_k(:,2)];
            tpsds=[tpsds,psds_k(:,1)];
            psds=[psds,psds_k(:,2)];

            PSs=[mean(tpss,2),mean(pss,2)];
              PSDs=[mean(tpsds,2),mean(psds,2)];
            IPSs=cumtrapz(PSDs(:,1),PSDs(:,2));
          

            %     ji=find(fx>(Fexc+20) & fx<(Fexc+30));
            ij=find(fx>(Fexc-2*Fexc_window) & fx<(Fexc-Fexc_window));
            ji=find(fx>(Fexc+Fexc_window) & fx<(Fexc+2*Fexc_window));
            pijs=polyfit(fx(ij),IPSs(ij),1);
            pjis=polyfit(fx(ji),IPSs(ji),1);
            Wdrive=sqrt(2*(polyval(pjis,Fexc)-polyval(pijs,Fexc)));
            Adrive1=Wdrive*k_exc;

            pij=polyfit(fx(ij),IPSx(ij),1);
            pji=polyfit(fx(ji),IPSx(ji),1);
            Wex=polyval(pji,Fexc)-polyval(pij,Fexc);

            figure, loglog(PSDs(:,1),IPSs,'b','linewidth',2)
            hold on, loglog(Fexc,polyval(pjis,Fexc),'*',Fexc,polyval(pijs,Fexc),'*')
            xlabel('f (Hz)'), ylabel('Excitation Spectrum $\int $P(f) ($V^2*s$)')
        end
    end
    
    PSx=[mean(tpsx,2),mean(psx,2)];
    PSDx=[mean(tpsdx,2),mean(psdx,2)];
    
    fx=PSDx(:,1);
    IPSx=cumtrapz(PSDx(:,1),PSDx(:,2));
    

    %     ij=find(fx>(Fexc-30) & fx<(Fexc-20));


    %     Adrive1=1e2*140*(64/16384)*0.02; disp(' SET Adrive, ln 146 !!!!')




    %figure, loglog(PSx(:,1),PSx(:,2),'b','linewidth',2)
    figure, loglog(PSx(:,1),PSx(:,2),'b','linewidth',2)
    xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')
    disp('Click on corner frequency (Hz) and steady-state amplitude ($V^2$)')
    [fc,Pss]=ginput(1);
    %    [fc,yo]=ginput(1);
    %    disp('Click on steady-state amplitude ')
    %    [xo,Pss]=ginput(1);

    Dfit=Pss*pi^2*(100^2+fc^2);
    gamma_bead=6*pi*viscosity*Rbead; %theoretical drag coefficient of bead
    %     fc=input('Corner frequency: ');
    Dtheor=kT/gamma_bead;
    Beta=sqrt(Dtheor/Dfit);

    %     Beta=130; % initial guesses
    %     fc=1e3;
    %     Dfit=15;


    z_flyv=0;%input('Use Flyvberg routine to fit thermal noise (0,1)? ');
    if z_flyv==1;
        % use calibration routine
        cd ./CPCpackage-v2
        crdat=[cposx,cposy];%zeros(Ncdat,2)];
        save calib_data.txt crdat -ascii -tabs
        %     start_fit
        disp('Load "calib_data.txt" with options: x=1, y=2, z=0, Fs=20, Nblock=100, 1 filter @ 10 kHz, fit 100<f<7000')
        run('./start_fit.m')
        zcont=input('Continue (enter after fitting power spectrum) (0/1)? ');
        cd ..
        fc=parametersX(1);
        Dfit=parametersX(2);
        %     Beta=sqrt(Pexc/(Aexc^2));

        Wth=Adrive1^2/(2*(1+(fc/Fexc)^2));
        Beta=sqrt(Wth/Wex);
    end
    %     calculate Px
    %         Pxexc=max(Px);
    %         iexc=find(Px==Pxexc);
    %         Pthexc=P_theor(scal_fit,parametersX,f(iexc),0);
    %         Wex=(Pxexc-Pthexc)*(1/ctime(end));
    %
    %         IPx=cumtrapz(f,Px);
    %         ii=find(f>(Fexc-1.5) & f<(Fexc-0.5));
    %         fi=find(f>(Fexc+0.5) & f<(Fexc+1.5));
    %         pii=polyfit(f(ii),IPx(ii),1);
    %         pfi=polyfit(f(fi),IPx(fi),1);
    %         Wex=polyval(pfi,Fexc)-polyval(pii,Fexc);
    %         %Wex=mean(IPx(fi))-mean(IPx(ii));
    %         Wth=(Aexc^2)/(2*(1+(fc/Fexc)^2)); %[nm^2*s] theoretical (see Tolic-Norrelykke, 2006, eq. 9)
    %
    %
    %     Beta=sqrt(Wth/Wex);
    %
    %
    %     falias=10e3;


    %     figure, loglog(PSx(:,1),PSx(:,2),'b',fx,Pth,'r','linewidth',2)
    %     xlabel('f (Hz)'), ylabel('Power Spectrum (V^2*s)')
    %     axis([10 1e4 min(Pth) max(Pth)])
    %
    %     figure, loglog(fx,IPSx,'b',fx,IPth,'r','linewidth',2)
    %     xlabel('f (Hz)'), ylabel('Cum. Power Spectrum (V^2)')
    %     axis([10 1e4 1e-5 max(IPth)+1e-2])
    %
    %     try increasing accuracy by fitting cumulative power spectrum
    disp('Now fit parameters to cumulative power spectrum --------')

    %     first to thermal response
    %Dfit=0.72;
    %fc=757;
    if f0==0,
        u10=[Dfit fc];% Qo*1e2 f0*1e2]; %initial guess
        par_range=[0.01 0.01];%[0.2*Dfit 0.1*fc];% 1e3 1e3];
        lb=[];%u10-par_range;%[0 500 50];
        ub=[];%u10+par_range;%[50 3e3 250];
    else
        u10=[Dfit fc Qo*1e6 f0];
        par_range=[0.2*Dfit 0.1*fc Qo f0];
        lb=[0,100,0,0];
        ub=[50,2000,1e-4,150];
    end
    disp([' - fit using ',num2str(numel(u10)),' parameters -'])
    options=optimset('TolFun',1e-12,'TolX',1e-10,'Display','on','MaxFunEvals',5e4,'MaxIter',1e4);%,...
    %'LargeScale','off','LevenbergMarquardt','on');%,'DiffMinChange',0.1);%,'DiffMinChange',5);%,'MaxFunEvals',1e3,'Display','iter');
    jfit1=find((fx>fmin & fx<fmax) & IPSx>0 & isnan(IPSx)~=1);
    jjf=find(fx(jfit1)<(Fexc-Fexc_window) | fx(jfit1)>(Fexc+Fexc_window));
    jfit1=jfit1(jjf);
    fdata1=fx(jfit1);%(Fexc+10):10:7e3;
    ydata1=cumtrapz(fx(jfit1),PSDx(jfit1,2));%interp1(fx,IPSx,fdata,'pchip')
    ydat10=ydata1;
    ydatr1=ydata1./ydat10;
    ydatr1(1)=1;
    %     ydatr=ydata;

    IPth1=ydat10.*excited_power_spectrum_thermal_v4(u10,fdata1,falias,Fexc,Adrive1,ydat10);
    figure, semilogx(fdata1,ydata1,'b',fdata1,IPth1,'r','linewidth',2)
    xlabel('f (Hz)'), ylabel('$P_{th}$ ($V^2*s$)')
    axis tight

    %options=optimset(options,'DiffMinChange',0.1,'TolX',1e-10,'TolFun',1e-12,'DiffMinChange',1);
    [u1,resnorm,residual]=lsqcurvefit(@excited_power_spectrum_thermal_v4,u10,fdata1,ydatr1,lb,ub,options,falias,Fexc,Adrive1,ydat10);

    IPth1=ydat10.*[1+residual];%excited_power_spectrum(u,fx(jfit),falias,Fexc,Adrive1,ydat0,f0);
    hold on, semilogx(fx(jfit1),IPth1,'g','linewidth',2)

    %%     now fit Beta
    if Fexc ~=0;
        %     mk=numel(cposx);
        %     %fsamp=100e3;
        %     Nfft=floor(Fs/10);%10e3;
        %     Noverlap=floor(Nfft/2);
        %     feps=Fexc/10;
        %     [P,F]=spectrum(cstagex,cposx,mk,Noverlap,hanning(Nfft),Fs,0.95);
        %     if_range=find(abs(F-Fexc)<feps);
        %     if_des=if_range(floor(median(find(P(if_range,5)==max(P(if_range,5))))));
        %     Fdes=F(if_des);
        %     Pxx=P(if_des,1);
        %     Pyy=P(if_des,2);
        %     H=P(if_des,4);
        %     C=P(if_des,5);
        %     Md1=abs(H(:));
        %     Wd1=2*pi*Fdes;
        %     Wc=2*pi*u1(2);
        %     Beta=(k_exc/Md1)*Wd1/sqrt(Wd1^2+Wc^2);
        %
        %
        jfit2=find((fx>fmin & fx<fmax) & IPSx>0 & isnan(IPSx)~=1);
        fdata2=fx(jfit2);%(Fexc+10):10:7e3;
        ydata2=cumtrapz(fx(jfit2),PSDx(jfit2,2));%interp1(fx,IPSx,fdata,'pchip')
        ydat20=ydata2;
        ydatr2=ydata2./ydat20;
        ydatr2(1)=1;

        u20=Beta;
        lb2=0;
        ub2=1e5;


        IPth2=ydat20.*excited_power_spectrum_v4(u20,fdata2,falias,Fexc,Adrive1,ydat20,u1);
        figure, semilogx(fdata2,ydata2,'b',fdata2,IPth2,'r','linewidth',2)
        xlabel('f (Hz)'), ylabel('$P_{th}+P_{ex}$ $(V^2*s)$')
        axis tight

        u2=Beta;
        %     Fexc=Fdes;
        options=optimset(options,'DiffMinChange',1,'TolX',1e-10,'TolFun',1e-12);
        [u2,resnorm,residual]=lsqcurvefit(@excited_power_spectrum_v4,u20,fdata2,ydatr2,lb2,ub2,options,falias,Fexc,Adrive1,ydat20,u1);

        IPth2=ydat20.*[1+residual];%excited_power_spectrum(u,fx(jfit),falias,Fexc,Adrive1,ydat0,f0);
        hold on, loglog(fx(jfit2),IPth2,'c--','linewidth',2)
        axis tight

        zup=1;%input('Update parameters based on this fit (0,1)? ');
        if zup==1,
            Lu1=numel(u1);
            if Lu1==2,
                Dfit=u1(1);
                fc=u1(2);
                Beta=u2;

                Pth_t=(Dfit./(pi^2.*(fx.^2 + fc^2))).*(falias^2./(falias^2 + fx.^2));
                df=abs(fx-Fexc);
                iexc1=find(df == min(df));
                Pth_r1=zeros(size(fx));
                Pth_r1(iexc1)=(Adrive1/Beta)^2./(2*(1+fc^2/Fexc^2));%Adrive^2./2
                Pth=Pth_t+Pth_r1;
                IPth=cumtrapz(fx,Pth);
            elseif Lu1==4,
                Dfit=u1(1);
                fc=u1(2);
                Beta=u2;
                Qo=1e-6*u1(3);
                f0=u1(4);

                Pth_t=((Dfit./(pi^2.*(fx.^2 + fc^2)))+(Qo*f0^2)./(f0^2+fx.^2)).*(falias^2./(falias^2 + fx.^2));
                df=abs(fx-Fexc);
                iexc1=find(df == min(df));
                Pth_r1=zeros(size(fx));
                Pth_r1(iexc1)=(Adrive1/Beta)^2./(2*(1+fc^2/Fexc^2));%Adrive^2./2
                Pth=Pth_t+Pth_r1;
                IPth=cumtrapz(fx,Pth);
            end
        end

        disp(['Corner frequency f_c = ', num2str(fc),' Hz'])
        disp(['Diffusion coefficient D = ',num2str(Dfit),' V = ',num2str(Beta^2*Dfit),' nm^2/s'])
        disp(['Beta = ',num2str(Beta),' nm/V'])
        k_trap=2*pi*fc*kT/(Dfit*Beta^2);
        gamma_fit=kT/(Dfit*Beta^2);
        disp(['Trap stiffness k_{trap} = ',num2str(k_trap),' pN/nm'])
        disp(['Drag coefficient \gamma = ',num2str(gamma_fit), ' pN*s/nm'])
        disp(['Beta * k_{trap} = ',num2str(Beta*k_trap),' pN/V'])


        disp('')
        PDsens_check=sqrt(Dtheor/Dfit);
        z_PDcheck=input(['Beta based on Dfit = ',num2str(PDsens_check),' nm/V, use this value instead (0/1)? ']);
        if z_PDcheck==1,
            Beta=PDsens_check;
            k_trap=2*pi*fc*kT/(Dfit*Beta^2);
            gamma_fit=kT/(Dfit*Beta^2);

            disp(['Trap stiffness k_{trap} = ',num2str(k_trap),' pN/nm'])
            disp(['Drag coefficient \gamma = ',num2str(gamma_fit), ' pN*s/nm'])
        end
    else
        
        Lu1=numel(u1);
        if Lu1==2,
            Dfit=u1(1);
            fc=u1(2);
            %Beta=u2;

            Pth=(Dfit./(pi^2.*(fx.^2 + fc^2))).*(falias^2./(falias^2 + fx.^2));
            IPth=cumtrapz(fx,Pth);
        elseif Lu1==4,
            Dfit=u1(1);
            fc=u1(2);
            %Beta=u2;
            Qo=1e-6*u1(3);
            f0=u1(4);

            Pth=((Dfit./(pi^2.*(fx.^2 + fc^2)))+(Qo*f0^2)./(f0^2+fx.^2)).*(falias^2./(falias^2 + fx.^2));
            IPth=cumtrapz(fx,Pth);
        end
    end
    PD_sensitivity=Beta;
    
        PDsens_check=sqrt(Dtheor/Dfit);
        disp(['Beta based on Dfit = ',num2str(PDsens_check),' nm/V']);
        Beta=PDsens_check;
        k_trap=2*pi*fc*kT/(Dfit*Beta^2);
        gamma_fit=kT/(Dfit*Beta^2);
        disp(['Corner frequency f_c = ', num2str(fc),' Hz'])
        disp(['Diffusion coefficient D = ',num2str(Dfit),' V = ',num2str(Beta^2*Dfit),' nm^2/s'])
        disp(['Trap stiffness k_{trap} = ',num2str(k_trap),' pN/nm'])
        disp(['Drag coefficient \gamma = ',num2str(gamma_fit), ' pN*s/nm'])


    figure, loglog(PSx(:,1),PSx(:,2),'b',fx,Pth,'r','linewidth',2)
    xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')
    %axis([10 1e4 min(Pth) max(Pth)])
    %     figure(4), hold on, loglog(Fexc,polyval(pji,Fexc),'b*')
    %     figure(4), hold on, loglog(Fexc,polyval(pij,Fexc),'b*')



    %     IPth=excited_power_spectrum([Dfit fc Beta],fx,falias,Fexc,Adrive1);
    %     figure, loglog(PSx(:,1),PSx(:,2),'b',fx,Pth,'r','linewidth',2)
    %     xlabel('f (Hz)'), ylabel('Power Spectrum ($V^2*s$)')
    %
    figure, semilogx(fx,IPSx,'b',fx,IPth,'r','linewidth',2)
    xlabel('f (Hz)'), ylabel('Cumulative Power Spectrum ($V^2$)')
    %axis([10 1e4 1e-5 max(IPth)+1e-2])


    % Debugging - mock excitation at 15 Hz
    %     ij=find(f>10 & f<14);
    %     ji=find(f>16 & f<20);
    %     pij=polyfit(f(ij),IPx(ij),1);
    %     pji=polyfit(f(ji),IPx(ji),1);
    %     W2ex=polyval(pji,f2)-polyval(pij,f2)
    %     W2th=(A2^2)/2;%(2*(1+fc^2/f2^2))
    %
    %     % check power spectrum fit
    %     Pth_t=(Dfit/pi^2)./(f.^2 + fc^2);
    %     Pth_r=zeros(size(f));
    %     df=abs(f-Fexc);
    %     iexc=find(df==min(df));
    %     Pth_r(iexc)=(Aexc/Beta)^2/(2*(1+(fc/Fexc)^2));
    %     Pth=Pth_t+Pth_r;
    %     IPth=cumtrapz(f,Pth);
    %
    %     figure, loglog(f,Px,'b',f,Pth,'r')
    %     figure, loglog(f,IPx,'b',f,IPth,'r')

    zsave=input('Save fit parameters (0/1)? ');
    if zsave==1,
        fsave=[datdir,'fit_parameters.txt'];
        fnum=str2num(calib_fname{kd}(13:15));
        fid=fopen(fsave,'a');
        fprintf(fid,'\n');
        fprintf(fid,'%g ',[fnum k_trap Beta k_trap*Beta fc Dfit gamma_fit f0 Qo Fs Rbead Fexc]);
        fclose(fid);
    end

end%

%% Calculations for the force traces starts here
if zcont==1, %analyze force trace data
    %load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cd('C:\Users\abdul\Desktop\All the data - LBC\optical trap data\160928-abdullah\test');
    fdat=textread(fullfile(datdir,force_fname{1}),'','headerlines',8);
    %fdat=fdat(1:500e3,:); disp('!!! truncated data set line 580 !!!')
    [Nfdat,Mfdat]=size(fdat);
    %ftime=[0:(Nfdat-1)]'./Fs_force; %[s]
    if cht==0,
        ftime=[0:(Nfdat-1)]'.*(1/Fs_force);
    else
        ftime=fdat(:,cht);%
    end
    
    if strcmp(zdim,'x')==1, %Mfdat<3,
        fpos=PD_sensitivity.*fdat(:,chypd);%1e3.*fdat(:,4); %[nm]
        fpos_off=PD_sensitivity.*fdat(:,chxpd);
    elseif strcmp(zdim,'y')==1,
        fpos=PD_sensitivity.*fdat(:,chxpd);
        fpos_off=PD_sensitivity.*fdat(:,chypd);
    elseif strcmp(zdim,'xy')==1,
        fxpos=PD_sensitivity.*fdat(:,chypd);
        fypos=PD_sensitivity.*fdat(:,chxpd);
        pxy = polyfit(fxpos,fypos,1);
        th = atan(pxy(1));
        Rmat = [cos(th), -sin(th); sin(th), cos(th)];
        fmt = [fxpos,fypos]*Rmat;
        fpos = fmt(:,1);
        fpos_off = fmt(:,2);

        figure, plot(ftime,fxpos,'b'), xlabel('Time (s)'), ylabel('Position (nm)')
        hold on, plot(ftime,fypos,'r'), legend('x','y')
        
        figure, plot(fxpos,fypos), xlabel('x (nm)'), ylabel('y (nm)')
        
        disp('Calculation based on Fx and Fy ...')
    end
    
    %
    %         % calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         gamma_bead=6*pi*viscosity*Rbead; %drag coefficient of bead
    %         %fc=input('Corner frequency: ');
    %         k_trap=2*pi*gamma_bead*fc;
    %         Dtheor=kT/gamma_bead;
    %         PDsens_check=sqrt(Dtheor/Dfit);
    %         z_PDcheck=input(['Update PD_sensitivity to ',num2str(PDsens_check),' (0/1)?']);
    %         if z_PDcheck==1, PD_sensitivity=PDsens_check; end

    %     force traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure, plot(ftime,fpos)
    xlabel('time (s)'), ylabel('distance (nm)')

    fpos0=0;%input('Enter zero-force position ("0" to choose from plot): ');
    if fpos0==0,
        disp('Choose reference window (zero force).')
        [x1,y1]=ginput(2);
        jref=find(ftime>x1(1) & ftime<x1(2));
        fpos0=mean(fpos(jref));
    end
    fpos=fpos-fpos0;

    ff=fpos.*k_trap;
    
    % Z-flip
    zflip=input('Flip plot (0/1)? ');
    if zflip==1%abs(min(ff))>max(ff),
        ff=-ff;
        fpos=-fpos;
    end

    %     median filter
    ftm=ftime;%medfilt1(ftime,Nmdn);
    [ffm,ffvar]=medvarfilt1(ff,Nmdn);
    ppm=ffm./k_trap;

    figure,
    plot(ftime,ff,'color','b');%[0.7 0.7 0.7])
    hold on, plot(ftm,ffm,'r','LineWidth',1)
    %axis([0 ftm(end) min(ffm) max(ffm)])
    xlabel('Time (s)'), ylabel('Force (pN)')
    %         figure, plot(ftm,ffvar,'Color',[0.7 0.7 0.7],'LineWidth',1), axis tight, grid on
    %         xlabel('time (s)'), ylabel('variance')

    %     find stall events
    z_stall_events=1;%input('Find stall events (0/1)? ');
    if z_stall_events==1,
        %disp('Set stall threshold.')
        %[xthreshold, Fthreshold]=ginput(1);
        Fthreshold=0.5; disp('Fthreshold = 0.5 pN')
        jstall=find(abs(ffm)>Fthreshold);
        jevent=find(diff(jstall)>10);
        jevent=[0;jevent;length(jstall)];
    
    for kevent=2:length(jevent);
        ji=jstall(jevent(kevent-1)+1);
        jf=jstall(jevent(kevent)); 
        stall_ti(kevent-1)=ftm(ji);
        stall_tf(kevent-1)=ftm(jf);
        stall_time(kevent-1)=ftm(jf)-ftm(ji);  % Plateau time
        ffm_segment=ffm(ji:jf); % Force trace segment 
        [stall_force_value,stall_force_index]=max(abs(ffm_segment)); % Index the maximum stall force value - could be negative or positive
        stall_force(kevent-1)=ffm_segment(stall_force_index); % Stall force
        jvel = ji:jf;   % Velocity segment
        pvel=polyfit(ftm(jvel),ppm(jvel),1);    % Velocity
        stall_vel(kevent-1)=pvel(1); % stall velocity
        jsnap=find(ftm>(ftm(jf)-0.05) & ftm<(ftm(jf)+0.05));
        psnap=polyfit(ftm(jsnap),ppm(jsnap),1);
        stall_snapvel(kevent-1)=psnap(1); %snap back velocity
    end

        %     apply conditions
        jkeep=find(stall_time > Cstall_time & abs(stall_vel) > Cstall_vel & abs(stall_snapvel) > Cstall_snapvel);
%         jkeep=find(stall_time > Cstall_time & abs(stall_vel) > Cstall_vel & abs(stall_snapvel) < Cstall_snapvel);
        
        for k=1:length(jkeep),
            hold on, plot([stall_ti(jkeep(k)),stall_tf(jkeep(k))],[Fthreshold,Fthreshold],'m-','LineWidth',2)
        end
        
        title(['Beta = ',num2str(PD_sensitivity),' nm/V, k_{trap} = ', num2str(k_trap),' pN/nm'])

        figure, hist(stall_force(jkeep)), xlabel('stall force (pN)'), ylabel('Count')

        stall_force_keep=stall_force(jkeep);
        zs=input('Save (0/1)? ');
        if zs==1,
            %fnsave1=fullfile('./trap-calibration-data/',[force_fname{1}(1:end-4),'in-vitro-force.mat']);
            fnsave1=fullfile(datdir,[force_fname{1}(1:end-4),'-force.mat']);
            save(fnsave1,'k_trap1','ffm','ftm','fpos','fpos_off','stall_force','jkeep','stall_time','stall_vel','stall_snapvel','stall_ti','stall_tf')
            fnsave=[force_fname{1}(1:end-4),'_stall_force.txt'];
            save([datdir,fnsave],'stall_force_keep','-ascii','-tabs')
        end
    end
    
    figure, plot(ftime,fpos,'b',ftm,ffm./k_trap,'r')
    xlabel('time (s)'), ylabel('position (nm)')
    
    figure
    [ax,h1,h2]=plotyy(ftime,fpos,ftime,ffm);
    ymax=max(fpos)+20;
    ymin=min(fpos)-20;
    set(ax(1),'YLim',[ymin ymax]);
    xlabel('Time (s)')
    set(get(ax(1),'YLabel'),'String','Position (nm)')
    %set(ax(1),'Ytick',[-400:8:400]), grid on
    set(ax(1),'Ytick',-400:100:400)
    set(ax(2),'YLim',k_trap.*[ymin ymax])
    set(ax(2),'Ytick',[-60:5:60])
    set(get(ax(2),'YLabel'),'String','Force (pN)')
    set(h2,'linewidth',2)
    set(gca,'Box','off')
    set(ax(1),'Ycolor','k')
    set(ax(2),'Ycolor','k')
    set(ax(1),'FontSize',30)
    set(ax(2),'FontSize',30)
    set(ax(1),'FontWeight','bold')
 
    
    %Emily trying some new stuff
    
%     figure
%     yyaxis left
%     plot(ftime, fpos);
%     ymax=max(fpos)+20;
%     ymin=min(fpos)-20;
%     ylim([ymin ymax]);
%     ylabel('Position (nm)','fontweight','bold','FontSize',24);
%     pbaspect([2 1 1]);
%     
%     yyaxis right
    figure
    ftime_smoothed=smooth(ftime, 10,'sgolay',1);
    ffm_smoothed=smooth(ffm,10,'sgolay',1);
    plot(ftime_smoothed, ffm_smoothed);
%     ylim([-25 25]);
    ylabel('Force (pN)','fontweight','bold','FontSize',24);
    xlabel('Time (s)','fontweight','bold','FontSize',24);
    pbaspect([2 1 1]);
    
    
    
%     set(get(ax(1),'YLabel'),'String','Position (nm)')
%     set(ax(2),'YLim',k_trap.*[ymin ymax])
%     set(ax(2),'Ytick',[-60:5:60])
%     set(get(ax(2),'YLabel'),'String','Force (pN)')
%     set(h2,'linewidth',2)
%     set(gca,'Box','off')
%     set(ax(1),'Ycolor','k')
%     set(ax(2),'Ycolor','k')
%     set(ax(1),'FontSize',30)
%     set(ax(2),'FontSize',30)
%     set(ax(1),'FontWeight','bold')
%  
    
    figure
    plot(ftime_smoothed,ffm_smoothed);
    ymax=25;
    ymin=-25;
    ylim([ymin ymax]);
    xlabel('Time (s)')
    hold on
    plot([0; 0], [-15; -10], '-k',  [0; 5], [-15; -15], '-k', 'LineWidth', 2)
    hold on
    plot([0; 60], [0;0], '-k')
    hold off
    text(-0.5,-13, '5 pN', 'HorizontalAlignment','right')
    text(2.5,-16, '5 s', 'HorizontalAlignment','center')
    set(gca, 'Visible', 'off')
    %End of Emily trying stuff
    figure, hist(stall_vel(jkeep)), xlabel('Velocity (nm/s)')
    
    

%     %     to plot all stall force data
%     dd=dir([datdir,'*stall_force.txt']);
%     %     jp=7:14;
%     FS=[];
%     for k=1:length(dd),
%         dk=textread([datdir,dd(k).name]);
%         FS=[FS,dk(1:(end-1))];
%     end
% 
%     [muhat,sigmahat]=normfit(FS);
%     ffit=(min(FS)-0.5):0.1:(max(FS)+0.5);
%     nfit=normpdf(ffit,muhat,sigmahat);
%     [Ndat,Fdat]=hist(FS);
%     figure, bar(Fdat,Ndat./sum(Ndat))
%     hold on, plot(ffit,0.7.*nfit,'r-','LineWidth',2)
%     xlabel('Stall Force (pN)'), ylabel('probability density')
%     %     disp('Choose location for text on plot.')
%     text(mean(ffit)+1.5, 0.3,{['N = ',num2str(sum(Ndat))],['\mu = ', num2str(muhat)],['\sigma = ',num2str(sigmahat)]})
% elseif zcont==2; %force feedback
%     fdat=textread([datdir,force_fname],'','headerlines',8);
%     [mf,nf]=size(fdat);
%     if nf==5,
%         zflip=input('Flip data (0/1)? ');
%         Nfdat=mf;
%         ftime=fdat(:,1);%[0:(Nfdat-1)]'./2e3; %[s]
%         fpos=PD_sensitivity.*fdat(:,chxpd);%1e3.*fdat(:,4); %[nm]
% 
%         if zflip~=1,
%             %             fpos=fpos;
%             fpid=-fdat(:,5);
%         else
%             fpos=-fpos;
%             fpid=fdat(:,5);
%         end
%         fforce=k_trap.*fpos;
%         figure, plot(ftime,fforce), xlabel('Time (s)'), ylabel('Force (pN)')
%         figure, hist(fforce,50), xlabel('Force (pN)'), ylabel('Counts')
% 
%         [fpm,fpvar]=medvarfilt1(fpid,Nmdn);
%         figure, plot(ftime,fpid,'Color',[0.4 0.4 0.4],'LineWidth',2)
%         hold on, plot(ftime,fpm,'Color',[0.8 0.8 0.8])
%         xlabel('Time (s)'), ylabel('PID output (V)')
%         title(['$<Force> = $',num2str(mean(fforce)),' $\pm$ ',num2str(std(fforce)),' pN'])
% 
%         k_aod=20e3/.5266; %[nm/V]
%         xpid=fpid.*k_aod;
%         figure, plot(ftime,xpid,'Color',[0.4 0.4 0.4])
%         xlabel('Time (s)'), ylabel('PID output (nm)')
%     end
end


%frequency response
% [P,F]=spectrum(vel,chdat,N_max,N_overlap,hanning(N_fft),params.fscan,0.95);
% %P=[Pxx Pyy Pxy Txy Cxy Pxxc Pyyc Pxyc]
% if_range=find(abs(F-f(i_f))<f_eps);
% if_des=if_range(floor(median(find(P(if_range,5)==max(P(if_range,5))))));
% Pxx(i_f,i_ch)=P(if_des,1);
% Pyy(i_f,i_ch)=P(if_des,2);
% H(i_f,i_ch)=P(if_des,4);
% C(i_f,i_ch)=P(if_des,5);

% Fs = 1000;
% t = 0:1/Fs:.3;
% x=cos(2*pi*t*200)+randn(size(t));
% window=33;
% noverlap=32;
% nfft=4097;
% h = spectrum.welch('Hann',window,100*noverlap/window);
% hpsd = psd(h,x,'NFFT',nfft,'Fs',Fs);
% Pw = hpsd.Data;
% Fw = hpsd.Frequencies;


%h=spectrum.welch('Hann',Nwindow,100*Noverlap/Nwindow);
%hpsd=psd(h,

% [P,F]=spectrum(cstagex,cposx,Nfft,0,hanning(Nfft),Fs,0.95);
% Pxx=P(:,1);
% Pyy=P(:,2);
% Pxy=P(:,3);
% Txy=P(:,4);
% Cxy=P(:,5);
%
% figure, subplot(211), semilogx(F,abs(Txy)), ylabel('Txy Mag.')
% subplot(212), semilogx(F,unwrap(angle(Txy))), xlabel('f, Hz'), ylabel('Txy Phase')



