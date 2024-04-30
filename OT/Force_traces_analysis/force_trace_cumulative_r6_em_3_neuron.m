% Penn              14 Feb 2012         Adam G. Hendricks
% aggregate force trace data from several experiments

%November 2021: Emily cleaning up some stuff that she added

%%
clear all, close all, clear java
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/histf_v1/');
% datdir='/Volumes/Emily_2022/Omar_OT/U2OS/Good_force_traces/Force_mat/';
datdir='/Volumes/Emily_2022/Neurons/Good_OT_motors/';
%flist=dir([datdir,'*_force.mat']); %LBCs in live cells
WT_color=[0.3 0.3 0.3]; %early
S421D_color=[0 0.6 0]; %early
% WT_color=[0.2902 0.5569  0.8863]; %late
% S421D_color=[0.4314    0.6863    0.4745];%late
for condition=1:2
    if condition==1
        flist=dir([datdir,'*.mat']);
    elseif condition ==2
        flist=dir([datdir,'*.mat']);
    end
 
numel(flist); %to print the number of files incorporated into the analysis

cmapa=colormap(winter(12));
cmapr=colormap(hot(12));
%%

stall_force=[];
stall_time=[];
stall_vel=[];
stall_snapvel=[];
time_ant=[];
time_ret=[];
max_force_ant=[];
max_force_ret=[];
Fthreshold =0.5; %[pN] force threshold to define and event - 1.5, Emily modified from 0.7 to 0.5, Emily added this from unbinding rate code

%Emily trying to incorporate binding rate
force_bins = [-15:0.5:-Fthreshold,Fthreshold:0.5:15];
C_unbinding = zeros(size(force_bins)); %number of unbinding events in each force bin
C_points = zeros(size(force_bins)); %number of data points in each force bin
binding_rate = [];

for kf=1:numel(flist)
    datf=load([datdir,flist(kf).name]);
    stall_force=[stall_force,datf.stall_force];
    stall_time=[stall_time,datf.stall_time];
    stall_vel=[stall_vel,datf.stall_vel];
    stall_snapvel=[stall_snapvel,datf.stall_snapvel];
    max_force_ant=[max_force_ant,max(datf.stall_force)];
    max_force_ret=[max_force_ret,min(datf.stall_force)];
    
    jkant=find(datf.stall_force>0);
    jkret=find(datf.stall_force<0);
    time_ant(kf)=sum(datf.stall_time(jkant));
    time_ret(kf)=sum(datf.stall_time(jkret));
    
    %Emily trying to add in the binding rate
    flat_part_C_points = hist(datf.stall_force,force_bins);
end

%apply stall criteria
%jkp=find(stall_time>0.1 & abs(stall_vel)>200);
%jkp=find(stall_time>0.1 & abs(stall_vel)>500);
%jkp=find(stall_time>1 & abs(stall_snapvel)>250 & abs(stall_vel)>0);% & sign(stall_snapvel)+sign(stall_vel)==0);%[0.25,1.0] & abs(stall_vel)>500);
% jkp=find(stall_time>1); %Original
jkp=find(stall_time>0.25); %Emily
%jkp=find(stall_time>0);

jkpa=find(stall_force(jkp)>0);
mean_force_ant=mean(stall_force(jkp(jkpa)));
std_force_ant=std(stall_force(jkp(jkpa)));

jkpr=find(stall_force(jkp)<0);
mean_force_ret=mean(stall_force(jkp(jkpr)));
std_force_ret=std(stall_force(jkp(jkpr)));

%Emily added to get true stall force values
mean_force_ant_2{condition}=mean_force_ant;
mean_force_ret_2{condition}=mean_force_ret;
std_force_ant_2{condition}=std_force_ant;
std_force_ret_2{condition}=std_force_ret;
stall_force_ant_2{condition}=stall_force(jkp(jkpa));
stall_force_ret_2{condition}=stall_force(jkp(jkpr));
stall_time_ant_2{condition}=stall_time(jkp(jkpa));
stall_time_ret_2{condition}=stall_time(jkp(jkpr));

% dfs=0.25;%[0.5]
% dfs=0.25; %Emily modifying based on Adam's suggestion
dfs=0.5; %Emily modifying based on Adam's suggestion-modified again to increase bin size
% fs=-15:dfs:15;%Original
fs=-25:dfs:25; %Emily
ns=hist(stall_force(jkp),fs);
figure %('Name',append(num2str(condition),'Force histogram 1'),'NumberTitle','off'), hold on, 
bar(fs,ns,1)
xlabel('Force (pN)'), ylabel('Fraction of Events')
%Emily trying some stuff
% figure('Name',append(num2str(condition),'Stall time force distribution'),'NumberTitle','off'), hold on,
% bar(stall_time(jkp), stall_force(jkp))
% xlabel('Force (pN)'), ylabel('Fraction of Events')

ja=find(fs>0);
jr=find(fs<0);

figure('Name',append(num2str(condition),'Force histogram 2'),'NumberTitle','off'), hold on, 
xlabel('Force (pN)'), ylabel('Fraction of Events'), legend('Anterograde','Retrograde')
hold on, hba=bar(fs(ja),ns(ja)./sum(ns(ja)),1) 
set(hba,'facecolor','w')%'edgecolor',[1 1 1],'linewidth',2)
hold on, hbr=bar(-fs(jr),-ns(jr)./sum(ns(jr)),1,'r')
set(hbr,'facecolor','w','edgecolor',0.25.*[1 1 1],'linewidth',2)
axis tight
nsa_max=max(ns(ja)./sum(ns(ja)));
nsr_max=max(ns(jr)./sum(ns(jr)));
publication_fig(0,0,1)


 %saving the variables
% force_fraction_ant{condition}=ns(ja)./(sum(ns(ja))+sum(ns(jr))); %Original
% force_fraction_ret{condition}=-ns(jr)./(sum(ns(ja))+sum(ns(jr))); %Original
force_fraction_ant{condition}=ns(ja)./(sum(ns(ja))+sum(ns(jr))); %Emily modified Mar 2022 to have total fraction rather than only fraction of anterograde
force_fraction_ret{condition}=-ns(jr)./(sum(ns(ja))+sum(ns(jr)));%Emily modified Mar 2022 to have total fraction rather than only fraction of retrograde
force_range_ant{condition}=fs(ja);
force_range_ret{condition}=-fs(jr);

% %fit histograms
% Nopt_max=8;%original
Nopt_max_ant=3; %Emily for late
% Nopt_max_ant=3; %Emily for early
fss_ant=0:0.01:20; %changed from 0:0.01:15
fss_ret=-20:0.01:0; %changed from -15:0.01:0
options = statset('Display','final','MaxIter',1000,'TolX',1e-8,'robust','on');

for kdist=1:Nopt_max_ant;%1:10; 
    objk_ant = gmdistribution.fit(stall_force(jkp(jkpa))',kdist,'Options',options,'Start','randSample');
    obj_force_ant{kdist}=objk_ant
    aic_force_ant(kdist)=objk_ant.AIC
    bic_force_ant(kdist)=objk_ant.BIC 
    pdfka=pdf(objk_ant,fss_ant');
    %hold on, plot(fss_ant,pdfka.*1.2.*(nsr_max/max(pdfkr)),'color',cmapa(kdist,:),'linewidth',2)
end

% Nopt_max_ret=2; %Emily for early
Nopt_max_ret=4; %Emily for late
for kdist=1:Nopt_max_ret;%1:10;
    objk_ret = gmdistribution.fit(stall_force(jkp(jkpr))',kdist,'Options',options,'Start','randSample');
    obj_force_ret{kdist}=objk_ret;
    aic_force_ret(kdist)=objk_ret.AIC;
    bic_force_ret(kdist)=objk_ret.BIC;
    pdfkr=pdf(objk_ret,fss_ret');
    %hold on, plot(-fss_ret,-pdfkr.*1.2.*(nsa_max/max(pdfkr)),'color',cmapr(kdist,:),'linewidth',2)
end
%[aic_llf_ret,jopt_ret]=min(aic_force_ret);
[bic_llf_ret,jopt_ret]=min(bic_force_ret);
pfit_stalls_ret=pdf(obj_force_ret{jopt_ret},fss_ret');
Cpr=1.2.*(max(ns(jr)./sum(ns(jr)))/max(pfit_stalls_ret))
pfit_stalls_ret_per_cond{condition}=pfit_stalls_ret.*Cpr; %Emily
hold on, plot(-fss_ret,-pfit_stalls_ret.*Cpr,'linewidth',2,'color',0.25.*[1 1 1])
% for kc=1:jopt_ret, %to plot each component separately
%     pfc=mvnpdf(fss_ret',obj_force_ret{jopt_ret}.mu(kc),obj_force_ret{jopt_ret}.Sigma(kc));
%     %pfc=pdf('norm',fss_ret',obj_force_ret{jopt_ret}.mu(kc),obj_force_ret{jopt_ret}.Sigma(kc));
%     hold on, plot(-fss_ret,-pfc.*Cpr*(obj_force_ret{jopt_ret}.PComponents(kc)),'k','linewidth',1)
% end

%[aic_llf_ant,jopt_ant]=min(aic_force_ant);
[bic_llf_ant,jopt_ant]=min(bic_force_ant);
pfit_stalls_ant=pdf(obj_force_ant{jopt_ant},fss_ant');
Cpa=1.2.*(max(ns(ja)./sum(ns(ja)))/max(pfit_stalls_ant));
pfit_stalls_ant_per_cond{condition}=pfit_stalls_ant.*Cpa; %Emily
hold on, plot(fss_ant,pfit_stalls_ant.*Cpa,'linewidth',2,'color',0.25.*[1 1 1])
% for kc=1:jopt_ant, %to plot each component separately
%     pfc=mvnpdf(fss_ant',obj_force_ant{jopt_ant}.mu(kc),obj_force_ant{jopt_ant}.Sigma(kc));
%     %pfc=pdf('norm',fss_ret',obj_force_ret{jopt_ret}.mu(kc),obj_force_ret{jopt_ret}.Sigma(kc));
%     hold on, plot(fss_ant,pfc.*Cpa*(obj_force_ant{jopt_ant}.PComponents(kc)),'k','linewidth',1)
% end

axis tight
set(gca,'xtick',[0:2:25])
pbaspect([2 1 1])

figure ('Name',append(num2str(condition),'Force histogram anterograde retrograde'),'NumberTitle','off'), hold on, 
plot(aic_force_ant,'b','linewidth',2)
hold on, plot(bic_force_ant,'b--','linewidth',2)
hold on, plot(aic_force_ret,'r','linewidth',2)
hold on, plot(bic_force_ret,'r--','linewidth',2)
ylabel('AIC, BIC')
xlabel('N components')

disp('Anterograde fits:')
obj_force_ant{jopt_ant}

disp('Retrograde fits:')
obj_force_ret{jopt_ret}

obj_ret_stalls = gmdistribution.fit(stall_force(jkp(jkpr))',2,'Options',options);
obj_ant_stalls = gmdistribution.fit(stall_force(jkp(jkpa))',2,'Options',options);
%k=1; objk=[obj_force_ant{k}.mu, obj_force_ant{k}.Sigma(:), obj_force_ret{k}.PComponents']; sortrows(objk)

figure ('Name',append(num2str(condition),'Mean force anterograde retrograde'),'NumberTitle','off'), hold on, 
hb00=errorbar(1,mean_force_ant,std_force_ant,'.','markersize',20,'linewidth',2,'color',[0 0.5 0])
hb01=errorbar(2,-mean_force_ret,std_force_ret,'.','markersize',20,'linewidth',2,'color',[0 0 0.5])
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'Anterograde','Retrograde'})
% errorbar_tick(hb00,0.5,'units') 
% errorbar(hb00,0.5,'units')
% errorbar_tick(hb01,0.5,'units')
% errorbar(hb01,0.5,'units')
ylabel('Mean Force (pN)')

df=5;
fsm=0:df:30;
nfsm_ant=hist(max_force_ant,fsm);
nfsm_ret=hist(-max_force_ret,fsm);
figure ('Name',append(num2str(condition),'Max force anterograde retrograde'),'NumberTitle','off'), hold on, 
bar(fsm-df/4,nfsm_ant,0.5,'b')
hold on, bar(fsm+df/4,nfsm_ret,0.5,'r')
xlabel('Max. Force per Trace (pN)'), ylabel('Traces')

%
F_time_ant=time_ant./(time_ant+time_ret);
F_time_ret=time_ret./(time_ant+time_ret);
mean_time_ant=mean(F_time_ant);
std_time_ant=std(F_time_ant);
sem_time_ant=std_time_ant/sqrt(numel(time_ant));

mean_time_ret=mean(F_time_ret);
std_time_ret=std(F_time_ret);
sem_time_ret=std_time_ret/sqrt(numel(time_ret));

figure ('Name',append(num2str(condition),'Force histogram ant ret 3'),'NumberTitle','off'), hold on, 
hb1=errorbar([1],[mean_time_ant],[std_time_ant],'.','markersize',20,'linewidth',2,'color','b')
hold on, hb2=errorbar(2,mean_time_ret,std_time_ret,'.','markersize',20,'linewidth',2,'color','r')
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'Anterograde','Retrograde'})
% errorbar_tick(hb1,0.5,'units')
% errorbar_tick(hb2,0.5,'units')
set(gca,'ylim',[0 1])
ylabel('Fraction of Time (s)')

jant=find(stall_force>0 & stall_time>0);
jret=find(stall_force<0 & stall_time>0);

ts=0:0.05:5;
figure ('Name',append(num2str(condition),'Fraction of events ant ret'),'NumberTitle','off'), hold on, 
nts_ant=hist(stall_time(jant),ts);
nts_ret=hist(stall_time(jret),ts);
figure, plot(ts,nts_ant./sum(nts_ant),'b-',ts,nts_ret./sum(nts_ret),'r-','linewidth',2)
xlabel('Time (s)'), ylabel('Fraction of Events'), legend('ant.','ret.')

figure ('Name',append(num2str(condition),'Force histogram stall force 4'),'NumberTitle','off'), hold on, 
 plot(stall_time(jant),stall_force(jant),'b.',stall_time(jret),-stall_force(jret),'r.')
xlabel('Time (s)'), ylabel('Max. Force (pN)')

figure ('Name',append(num2str(condition),'Force histogram stall force 5'),'NumberTitle','off'), hold on, 
scatter_ant=scatter(stall_time(jant),stall_force(jant),'MarkerFaceColor',[0.085 0.77 0],'MarkerEdgeColor','none'),hold on,
alpha(scatter_ant,0.2)
scatter_ret=scatter(stall_time(jret),-stall_force(jret),'b','filled')
alpha(scatter_ret,0.2)
xlabel('Time (s)'), ylabel('Max. Force (pN)')
xlim([0 5])
publication_fig(0,0,1)


%%
%Emily assigning variables for plotting later
stall_time_ant{condition}=stall_time(jant);
stall_time_ret{condition}=-stall_time(jret);
stall_force_ant{condition}=stall_force(jant);
stall_force_ret{condition}=stall_force(jret);

%%
%bootstrap to find confidence intervals
% [nbic_ant,naic_ant]=bootstrap_multicomp_fits(250,Nopt_max,stall_force(jkp(jkpa)));
% [nbic_ret,naic_ret]=bootstrap_multicomp_fits(250,Nopt_max,stall_force(jkp(jkpr)));
% 
% xnc=1:Nopt_max;
% 
% nnca=hist(nbic_ant,xnc);
% figure, hba2=bar(xnc,nnca./sum(nnca))
% set(hba2,'facecolor',0.25.*[1 1 1],'edgecolor',[1 1 1],'linewidth',2)
% xlabel({'Optimal Number of Components','Anterograde Forces'})
% 
% nncr=hist(nbic_ret,xnc);
% figure, hbr2=bar(xnc,nncr./sum(nncr))
% set(hbr2,'facecolor',[1.0 1.0 1.0],'edgecolor',0.25.*[1 1 1],'linewidth',2)
% xlabel({'Optimal Number of Components','Retrograde Forces'})
%%

% test that anterograde and retrograde force event histograms are distinct
% using kolomogrov-smirnov test

[h_ks,p_ks,k_ks]=kstest2(stall_force(jkp(jkpr)),stall_force(jkp(jkpa)))


%%
% STEPS %%%%%%%%%%%%%%

%flist_steps=dir([datdir,'*_ant_steps.mat']);
%flist_steps=dir([datdir,'*_ret_steps.mat']);
for js=1:2,
    if js==1,
        flist_steps=dir([datdir,'*_ant_steps.mat']);
        %flist_steps=dir([datdir,'*_ant_steps_sub5pN.mat']);
    elseif js==2,
        flist_steps=dir([datdir,'*_ret_steps.mat']);
        %flist_steps=dir([datdir,'*_ret_steps_sub5pN.mat']);
    end
    step_size=[];
    step_dwell=[];
    for kf=1:numel(flist_steps),
        datf=load([datdir,flist_steps(kf).name]);
        step_size=[step_size,datf.steps.steps(5,:)];
        dwellk=(1/2e3).*[datf.steps.steps(3,1),diff(datf.steps.steps(3,:))];
        step_dwell=[step_dwell,dwellk];
    end

    dx=[-24:2:24];
    ndx=hist(step_size,dx);

    %fit to gaussian
    % jda=find(dx>0);
    % jdr=find(dx<0);
    % jfa=find(step_size>0);
    % jfr=find(step_size<0);
    % [mu_a,s_a]=normfit(step_size(jfa));
    % [mu_r,s_r]=normfit(step_size(jfr));
    for case_fits=1:4,
        switch case_fits
            case 1
                obj_steps0(1).mu=[-8; 8];
                obj_steps0(1).Sigma(:,:,1)=8;
                obj_steps0(1).Sigma(:,:,2)=8;
                obj_steps0(1).PComponents=[0.5 0.5];
                ks=2;
            case 2
                obj_steps0(2).mu=[-8; -4; 8];
                obj_steps0(2).Sigma(:,:,1)=8;
                obj_steps0(2).Sigma(:,:,2)=8;
                obj_steps0(2).Sigma(:,:,3)=8;
                obj_steps0(2).PComponents=[0.4 0.2 0.4];
                ks=3;
            case 3
                obj_steps0(3).mu=[-8; 4; 8];
                obj_steps0(3).Sigma(:,:,1)=8;
                obj_steps0(3).Sigma(:,:,2)=8;
                obj_steps0(3).Sigma(:,:,3)=8;
                obj_steps0(3).PComponents=[0.4 0.2 0.4];
                ks=3;
            case 4
                obj_steps0(4).mu=[-8; -4; 4; 8];
                obj_steps0(4).Sigma(:,:,1)=8;
                obj_steps0(4).Sigma(:,:,2)=8;
                obj_steps0(4).Sigma(:,:,3)=8;
                obj_steps0(4).Sigma(:,:,4)=8;
                obj_steps0(4).PComponents=[0.4 0.1 0.1 0.4];
                ks=4;
        end 
        
%         objc = gmdistribution.fit(step_size',ks,'Options',options,'Start',obj_steps0(case_fits));
%         obj_steps{case_fits}=objc;
%         aic_steps(case_fits)=objc.AIC;
%         bic_steps(case_fits)=objc.BIC;
%     end
%     jbic=find(bic_steps==min(bic_steps));
% 
%     obj_steps{jbic}.mu
%     obj_steps{jbic}.Sigma
%     obj_steps{jbic}.PComponents
%     u0=[-8 3 17 8 3 48]; %ant
%     u0=[-8 3 98 8 3 35]; %ret

%     figure('Name','Step size 1','NumberTitle','off'), hold on, 
%     hbs=bar(dx,ndx,1), xlabel('Step Size (nm)'), ylabel('Counts')
%     set(hbs,'facecolor',0.25.*[1.0 1.0 1.0],'edgecolor',1.*[1 1 1],'linewidth',2)
%     ddx=dx(1):0.1:dx(end);
%     pdf_fit=pdf(obj_steps{jbic},ddx');
% %     hold on, plot(ddx,normpdf(ddx,mu_a,s_a).*(max(ndx(jda)))./(1/(2.5*s_a)),'c-','linewidth',2)
% %     hold on, plot(ddx,normpdf(ddx,mu_r,s_r).*(max(ndx(jdr)))./(1/(2.5*s_r)),'c-','linewidth',2)
%     hold on, plot(ddx,pdf_fit.*(max(ndx)./max(pdf_fit)),'color',0.25.*[1 1 1],'linewidth',2)
%     set(gca,'xtick',[-24:4:24])
%     
%      
%     xdwell=0:0.025:2.5;
%     ndwell=hist(step_dwell,xdwell);
%     
%     [phat,pci]=gamfit(step_dwell)
%     pdwell_fit=gampdf(xdwell,phat(1),phat(2));
%    
%     figure('Name','Dwell time','NumberTitle','off'), hold on, 
%     hbd=bar(xdwell,ndwell./sum(ndwell),1)
%     set(hbd,'facecolor',0.25.*[1 1 1],'edgecolor',1.*[1 1 1],'linewidth',0.5)
%     xlabel('Dwell Time (s)'), ylabel('Probability')
%     
%     hold on, plot(xdwell,pdwell_fit./sum(pdwell_fit),'k-')
%     set(gca,'xlim',[0 1])

end
%%

% % function res=N_normfit(u,x,n)
% %     %u=[mu1, sigma1, a1, mu2, sigma2, a2, ...]
% %     Nmodes=numel(u)/3;
% %     nfit=zeros(size(x));
% %     for kn=1:Nmodes
% %         uk=(kn-1)*3+[1:3];
% %         nk=normpdf(x,uk(1),uk(2));
% %         nfit=nfit+uk(3).*nk;
% %     end
% %     res=sum((n-nfit).^2);
end

%decide optimal number of distributions - see aicbic.m

end

% mean_ant_WT=mean(stall_force_ant{1});
% mean_ant_S421D=mean(stall_force_ant{2});
% mean_ret_WT=-mean(stall_force_ret{1});
% mean_ret_S421D=-mean(stall_force_ret{2});
% 
% sem_ant_WT=std(stall_force_ant{1})/(sqrt(numel(stall_force_ant{1})));
% sem_ant_S421D=std(stall_force_ant{2})/(sqrt(numel(stall_force_ant{2})));
% sem_ret_WT=std(stall_force_ret{1})/(sqrt(numel(stall_force_ret{1})));
% sem_ret_S421D=std(stall_force_ret{2})/(sqrt(numel(stall_force_ret{2})));

mean_ant_WT=mean_force_ant_2{1};
mean_ant_S421D=mean_force_ant_2{2};
mean_ret_WT=-mean_force_ret_2{1};
mean_ret_S421D=-mean_force_ret_2{2};

sem_ant_WT=std_force_ant_2{1};
sem_ant_S421D=std_force_ant_2{2};
sem_ret_WT=std_force_ret_2{1};
sem_ret_S421D=std_force_ret_2{2};

%Force histograms for all stalls above the threshold 0.25s stall force,
%binned by 0.5pN bins 
figure('Name',append(num2str(condition),'Force histogram both conditions'),'NumberTitle','off'), hold on,
xlabel('Force (pN)'), ylabel('Fraction of Events'),
hold on, 
b1=bar(force_range_ant{2},force_fraction_ant{2},1,'facecolor',S421D_color,'edgecolor',S421D_color,'linewidth',0.5);
b1.FaceAlpha=0.3;
hold on, 
b2=bar(force_range_ret{2},force_fraction_ret{2},1,'facecolor',S421D_color,'edgecolor',S421D_color,'linewidth',0.5);
b2.FaceAlpha=0.3;
hold on,
b3=bar(force_range_ant{1},force_fraction_ant{1},1,'facecolor',WT_color,'edgecolor',WT_color,'linewidth',0.5);
b3.FaceAlpha=0.3;
hold on, 
b4=bar(force_range_ret{1},force_fraction_ret{1},1,'facecolor',WT_color,'edgecolor',WT_color,'linewidth',0.5);
b4.FaceAlpha=0.3;
hold on, scatter(mean_ant_WT,0.085,100,'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
hold on, plot([mean_ant_WT-sem_ant_WT; mean_ant_WT+sem_ant_WT], [0.085; 0.085],'Color', WT_color);
hold on, scatter(mean_ret_WT,-0.085,100,'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
hold on, plot([mean_ret_WT-sem_ret_WT; mean_ret_WT+sem_ret_WT], [-0.085; -0.085],'Color', WT_color);
hold on, scatter(mean_ant_S421D,0.09,100,'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on, plot([mean_ant_S421D-sem_ant_S421D; mean_ant_S421D+sem_ant_S421D], [0.09; 0.09], 'Color', S421D_color);
hold on, scatter(mean_ret_S421D,-0.09,100,'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on, plot([mean_ret_S421D-sem_ret_S421D; mean_ret_S421D+sem_ret_S421D], [-0.09; -0.09],'Color', S421D_color);
xlim([-1.5 25])
publication_fig(0,0,1)
% axis tight

%Scatter plot all stalls, not filtered for stall >0.25s
y=zeros(6);
figure('Name',append(num2str(condition),'Force histogram stall force scatter both conditions'),'NumberTitle','off'); 
hold on 
scatter_ant_WT=scatter(stall_time_ant{1},stall_force_ant{1},'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ant_WT,0.3)
scatter_ant_S421D=scatter(stall_time_ant{2},stall_force_ant{2},'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ant_S421D,0.3)
scatter_ret_WT=scatter(-stall_time_ret{1},stall_force_ret{1},'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
alpha(scatter_ret_WT,0.3)
scatter_ret_S421D=scatter(-stall_time_ret{2},stall_force_ret{2},'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ret_S421D,0.3);
plot(0:5,y,'k'), hold on,
% plot(repelem(0.1,50),-49:2:49,'b.'),hold on,
% plot(repelem(0.25,50),-49:2:49,'.','MarkerEdgeColor',[0.87 0.49 0])
xlabel('Time (s)'), ylabel('Max. Force (pN)')
xlim([0 5])
ylim([-25 25])
% yticks([-50 -25 0 25 50])
publication_fig(0,0,1)
% pbaspect([2 1 1])

%Scatter plot for all stalls above threshold of 0.25s
y=zeros(6);
figure('Name',append(num2str(condition),'Force histogram stall force scatter both conditions v2'),'NumberTitle','off'); 
hold on 
scatter_ant_WT=scatter(stall_time_ant_2{1},stall_force_ant_2{1},'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ant_WT,0.3)
scatter_ant_S421D=scatter(stall_time_ant_2{2},stall_force_ant_2{2},'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ant_S421D,0.3)
scatter_ret_WT=scatter(stall_time_ret_2{1},stall_force_ret_2{1},'MarkerFaceColor',WT_color,'MarkerEdgeColor','none');
alpha(scatter_ret_WT,0.3)
scatter_ret_S421D=scatter(stall_time_ret_2{2},stall_force_ret_2{2},'MarkerFaceColor',S421D_color,'MarkerEdgeColor','none');
hold on
alpha(scatter_ret_S421D,0.3);
plot(0:5,y,'k'), hold on,
% plot(repelem(0.1,50),-49:2:49,'b.'),hold on,
% plot(repelem(0.25,50),-49:2:49,'.','MarkerEdgeColor',[0.87 0.49 0])
xlabel('Time (s)'), ylabel('Max. Force (pN)')
xlim([0 5])
ylim([-25 25])
% yticks([-50 -25 0 25 50])
publication_fig(0,0,1)
% pbaspect([2 1 1])

%Emily trying the bootstrapping to get a sortof sliding average with
%confidence intervals
% force_time_WT_ant=[stall_force_ant{1}; stall_time_ant{1}].';
% force_time_S421D_ant=[stall_force_ant{2}; stall_time_ant{2}].';
% force_time_WT_ret=[abs(stall_force_ret{1}); abs(stall_time_ret{1})].';
% force_time_S421D_ret=[abs(stall_force_ret{2}); abs(stall_time_ret{2})].';
force_time_WT_ant=[stall_time_ant_2{1};stall_force_ant_2{1}].';
force_time_S421D_ant=[stall_time_ant_2{2}; stall_force_ant_2{2}].';
force_time_WT_ret=[abs(stall_time_ret_2{1}); abs(stall_force_ret_2{1})].';
force_time_S421D_ret=[abs(stall_time_ret_2{2}); abs(stall_force_ret_2{2})].';



Loic_bootstrap_code_04092019_OT_em(force_time_WT_ant, 1000,0.05,3,5)  %Emily can't figure out how to make it work for plots with the same size bins
Loic_bootstrap_code_04092019_OT_em(force_time_S421D_ant, 1000,0.05,3,5)
Loic_bootstrap_code_04092019_OT_em(force_time_WT_ret, 1000,0.05,3,5)
Loic_bootstrap_code_04092019_OT_em(force_time_S421D_ret, 1000,0.05,3,5)


% % 
% fileprefix=append('20220404_late_stall_dfs_0pt5pN_0pt5s_fit_maxant_3_maxret_4_',num2str(condition))
% % % 
% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end

%do kstest here
% [kstest_ant,kstest_ant_p]=kstest2(anterograde_forces{1},anterograde_forces{2}) %anterograde forces
% 
% [kstest_ret,kstest_ret_p]=kstest2(retrograde_forces{1},retrograde_forces{2})
% 
% [kstest_ant_ret_WT,kstest_ant_ret_WT_p]=kstest2(anterograde_forces{1},retrograde_forces{1})
% 
% [kstest_ant_ret_S421D,kstest_ant_ret_S421D_p]=kstest2(anterograde_forces{2},retrograde_forces{2})

[kstest_ant,kstest_ant_p]=kstest2(force_fraction_ant{1},force_fraction_ant{2}) %anterograde forces

[kstest_ret,kstest_ret_p]=kstest2(abs(force_fraction_ret{1}),abs(force_fraction_ret{2}))

[kstest_ant_ret_WT,kstest_ant_ret_WT_p]=kstest2(force_fraction_ant{1},abs(force_fraction_ret{1}))

[kstest_ant_ret_S421D,kstest_ant_ret_S421D_p]=kstest2(force_fraction_ant{2},abs(force_fraction_ret{2}))