% This code computes the dynamics of the model.

%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');


%% Parameter Vector from Static Model
% load('Estimate_static');
load('GE_Results')

muT=bestParams(1);
sigmaT=bestParams(2);
tau=bestParams(3);
theta=bestParams(4);
F=bestParams(5);

sigma = 5;
% muT = 0.1345;
% sigmaT = 1.4128;
% tau = 1.3444;
% kappa =
% theta = 4.3136;
% sigma = 5;
% 
% Y = 126.7572;
% YF = 194.4945;
% LF = 173.4218;
% sigma=5;
% theta=kappa*(sigma-1);
% f=f*4.93*.43e-5;
% F=f/sigma;

% F = 1.0202e-5;
% New parameters for dynamic model
nu = 0.045;
mu = -theta*nu^2/2;
rho = 0;



%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Number of shadow firms
MF_large=10000;
MF_small=1400;
MH_large=5000;
MH_small=700;

% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Scale of model (i.e. number of sectors)
scale = 1;

% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4*scale;                                         
S_init = length(cdshares_init);                         % Number of sectors in data
S = S_init*S_multiple;                                  % Number of sectors used (see footnote 56)

% Assign CD-shares across sectors 
ALPHA = cdshares_init;      

for iloop = 1:S_multiple-1
    ALPHA = [ALPHA;cdshares_init(randperm(S_init))];
end
ALPHA = ALPHA/S_multiple;

% Split sectors into large and small based on CD-share.
split_param = 1.25;
small = (ALPHA<split_param/S);                            
Nsmall = sum(small);                           
Nlarge = S-Nsmall;

% Make random draws for loop
rng(aseed);                                               % Reset random number generator for consistency with old code
rtdraws = randn(1,S);

UH0S = exprnd(1,MH_small,Nsmall);                         % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1
UF0S = exprnd(1,MF_small,Nsmall);                        
UHS = cumsum(UH0S);                                       % Cumulate to get U of each home firm in each sector (see footnote 57)
UFS = cumsum(UF0S);

UH0L = exprnd(1,MH_large,Nlarge);
UF0L = exprnd(1,MF_large,Nlarge);
UHL = cumsum(UH0L);
UFL = cumsum(UF0L);

% Labor Normalization
L0 = 100;

% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z 
RT=exp(muT+sigmaT*rtdraws);
RTS=RT(small);
RTL=RT(~small);

% Draw productivities phi (step 2 of estimation procedure)
ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
ZFS=UFS.^(-1/theta);
ZHL=(UHL./(repmat(RTL,MH_large,1))).^(-1/theta);                        
ZFL=UFL.^(-1/theta);

%% Simulate Dynamics

% Length of dynamic simulation

% ZHS_Save = [];
rng(aseed);

% Years that are actually being recorded
RECORD = [1:11];
R_length = length(RECORD);
T = RECORD(end);

% Determine productivity draws
% su = sqrt(rho*nu^2);
% sv = sqrt(nu^2*(1-rho));

counter = 1;

% Compute Barrier: Minimum over all draws
min_HS=min(ZHS); %%% barrier= minimum over all draws
min_FS=min(ZFS);
min_HL=min(ZHL);
min_FL=min(ZFL);

[K,KF,PHI,PHIF,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,PHIFVEC,mom,varphi_bar,X,D,DSHM_small,DSHM_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
disp('Loop 0 is finished')
% save('Data/GE_Results','bestParams','Y','YF','LF','varphi_BB')
[Momarray]=Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF); 

for t=1:T
%     uz_HS = repmat(randn(1,S),ZHS_length,1);
%     uz_HL = repmat(randn(1,S),ZHL_length,1);
%     uz_FS = repmat(randn(1,S),ZFS_length,1);
%     uz_FL = repmat(randn(1,S),ZFL_length,1);
% 
%     vz_HS = randn(ZHS_length,S);
%     vz_HL = randn(ZHL_length,S);
%     vz_FS = randn(ZFS_length,S);
%     vz_FL = randn(ZFL_length,S);
%     
%     eps_HS = su*uz_HS+sv*vz_HS;
%     eps_HL = su*uz_HL+sv*vz_HL;
%     eps_FS = su*uz_FS+sv*vz_FS;
%     eps_FL = su*uz_FL+sv*vz_FL;
    rng(t);
    ZHS = exp( repmat(log(min_HS),MH_small,1) + abs( log(ZHS) - repmat(log(min_HS),MH_small,1) + mu + nu*randn(Nsmall,MH_small)' ));
    ZFS = exp( repmat(log(min_FS),MF_small,1) + abs( log(ZFS) - repmat(log(min_FS),MF_small,1) + mu + nu*randn(Nsmall,MF_small)' ));
    
    ZHL = exp( repmat(log(min_HL),MH_large,1) + abs( log(ZHL) - repmat(log(min_HL),MH_large,1) + mu + nu*randn(Nlarge,MH_large)' ));
    ZFL = exp( repmat(log(min_FL),MF_large,1) + abs( log(ZFL) - repmat(log(min_FL),MF_large,1) + mu + nu*randn(Nlarge,MF_large)' ));
    
    
    
%     ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+nu*randn(length(ZHS(~inactive_HS)),1)));
% %     ZHS_Save =[ZHS_Save,ZHS(:)];
%     ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+nu*randn(length(ZHL(~inactive_HL)),1)));
%     ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+nu*randn(length(ZFS(~inactive_FS)),1)));
%     ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+nu*randn(length(ZFL(~inactive_FL)),1)));
%     
%     ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+eps_HS(~inactive_HS)));
% %     ZHS_Save =[ZHS_Save,ZHS(:)];
%     ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+eps_HL(~inactive_HL)));
%     ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+eps_FS(~inactive_FS)));
%     ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+eps_FL(~inactive_FL)));
%     
    % If the year is part of RECORD, record PE results
    if any(RECORD==t)
        counter = counter+1;
        [K,KF,PHI,PHIF,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,PHIFVEC,mom,varphi_bar,X,D,DSHM_small,DSHM_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);        LAMBDAFVEC_t(counter,:) = LAMBDAFVEC;
        PHIFVEC_t(counter,:) = PHIFVEC;
        save_share{t} = DSHM_small';
        save_share_big{t} = DSHM_large';
        disp(['Loop ',num2str(t),' is finished']);
    end
    
end

% %% MEAN REVERSION
% 
% % GAMMAF
% GAMMAFVEC_t = LAMBDAFVEC_t-PHIFVEC_t;
% 
% % Delta LAMBDAF
% DELTA_LAMBDAFVEC_t = LAMBDAFVEC_t-repmat(LAMBDAFVEC_t(1,:),R_length+1,1);
% 
% % Prepare plotting variables
% ptiles = 10;
% LPCT = prctile(GAMMAFVEC_t(1,:),[1/ptiles:1/ptiles:1-1/ptiles]*100)';
% 
% mean_20 = zeros(ptiles,1);         
% mean_50 = zeros(ptiles,1);
% 
% min_20 = zeros(ptiles,1);
% min_50 = zeros(ptiles,1);
% 
% max_20 = zeros(ptiles,1);
% max_50 = zeros(ptiles,1);
% 
% DELTA_20 = DELTA_LAMBDAFVEC_t(2,:);
% DELTA_50 = DELTA_LAMBDAFVEC_t(3,:);
% 
% mean_20(1) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)));
% mean_50(1) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)));
% min_20(1) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)),10);
% min_50(1) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)),10);
% max_20(1) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)),90);
% max_50(1) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)),90);
% for i=2:ptiles-1
%    mean_20(i) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
%    mean_50(i) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
%    min_20(i) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),10); 
%    min_50(i) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),10); 
%    max_20(i) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),90); 
%    max_50(i) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),90); 
% end
% mean_20(ptiles) = mean(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
% mean_50(ptiles) = mean(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
% min_20(ptiles) = prctile(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),10);
% min_50(ptiles) = prctile(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),10);
% max_20(ptiles) = prctile(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),90);
% max_50(ptiles) = prctile(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),90);
% 
% 
% % Plot mean reversion graph
% figure(5)
% bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
% hold on
% bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
% lgd = legend('50 years','20 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.08, .065])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion.png')
% 
% figure(6)
% hold on
% bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
% errorbar(1:10,mean_20,min_20,max_20,'.','LineWidth',2,'Color',[0, 0.4470, 0.7410])
% lgd = legend('20 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.2, .2])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion_20.png')
% 
% figure(7)
% hold on
% bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
% errorbar(1:10,mean_50,min_50,max_50,'.','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
% lgd = legend('50 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.3, .3])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion_50.png')


% %% Computing SR and LR persistence (Table 4)
% 
% keepstd = zeros(S,T-1)*NaN;
% keepcorr = zeros(S,1)*NaN;
% % First the small sectors
% for z=1:Nsmall
%     
%     % Standard deviations across firms within time-sector
%     for tt=1:T-1
%         a2 = save_share_small{tt+1}(:,z);
%         a1 = save_share_small{tt}(:,z);
%         idx = (a2>0).*(a1>0);               % Only if firm exists in both periods
%         keepstd(z,tt) = std(a2(idx>0)-a1(idx>0));
%     end
%     
%     % Correlation over firms
%     a0 = save_share_small{1}(:,z);
%     aT = save_share_small{T}(:,z);
%     idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
%     if sum(idx)>1
%         keepcorr(z) = corr(aT(idx>0),a0(idx>0));
%     end
% end
% 
% % Now do the same for large sectors
% 
% for z=1:Nlarge
%     
%     for tt=1:T-1
%         a2 = save_share_large{tt+1}(:,z);
%         a1 = save_share_large{tt}(:,z);
%         idx = (a2>0).*(a1>0);
%         keepstd(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0));
%     end
%     
%     a0 = save_share_large{1}(:,z);
%     aT = save_share_large{T}(:,z);
%     idx = (a0>0).*(aT>0);
%     if sum(idx)>1
%         keepcorr(Nsmall+z) = corr(aT(idx>0),a0(idx>0));
%     end
%          
% end
% 
% LR_persistence = nanmedian(keepcorr);
% SR_persistence = nanmedian(nanmedian(keepstd,2));
% 
% DynMom = [SR_persistence,LR_persistence];

Nloop=T;

for i=1:Nsmall
    
    
    for tt=1:Nloop-1
        a2=save_share{tt+1}(i,:);
        a1=save_share{tt}(i,:);
        idx=(a2>0).*(a1>0);
        keepstd(i,tt)=std(a2(idx>0)-a1(idx>0));
    end
    
    %
    a0=save_share{1}(i,:);
    
    a8=save_share{Nloop}(i,:);
    idx=(a8>0).*(a0>0);
    if sum(idx)>1
        keepcorr(i)=corr(a8(idx>0)',a0(idx>0)');
    end
end

for i=1:S-Nsmall
    
    
    for tt=1:Nloop-1
        a2=save_share_big{tt+1}(i,:);
        a1=save_share_big{tt}(i,:);
        idx=(a2>0).*(a1>0);
        keepstd(Nsmall+i,tt)=std(a2(idx>0)-a1(idx>0));
        % standard deviation of delta rate of market share between two
        % years, by sector
    end
    %
    a0=save_share_big{1}(i,:);
    
    a8=save_share_big{Nloop}(i,:);
    idx=(a8>0).*(a0>0);
    if sum(idx)>1
        keepcorr(Nsmall+i)=corr(a8(idx>0)',a0(idx>0)');
    end
end





%% Compare with data moments

%%% Note = measure  MEDIAN across sectors measures
%%% play with nu until match the following:


longrun=median(keepcorr);
%%% French data/ stata goal 0.8678
shortrun=nanmedian(nanmedian(keepstd,2));