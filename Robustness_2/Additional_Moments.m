% This code computes additional moments given a vector of parameters. 
% It computes them in PE with 10K+ sectors.

%% Housekeeping
clear;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');
addpath('Results');

%% Load Parameters
load('Estimate_seed1_grid2') 

% Extract parameters
muT = bestParams(1,1);
sigmaT = bestParams(1,2);
tau = bestParams(1,3);
kappa = bestParams(1,4);
f = bestParams(1,5);

% Check that there is just one optimal parameter
if size(bestParams,1)~= 1
    warning('There is more than one optimal parameter combination.')
end

% Fix sigma
sigma = 5;

% Compute derivative parameters
theta = (sigma-1)*kappa;
F = f*4.93*.43*10^(-5)/5;

%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Relative size of economies
k = 1.75;

% Parameter governing the number of firms in each sector
M = 350; 

% Productivities for added firms
epsilon = 1e-10;

% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4;                                         
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

% Labor Normalization
L0 = 100;

% Get number of draws
MH = round(M*S*ALPHA);          
MF = round(k*M*S*ALPHA);        

mu_H = muT+sigmaT*rtdraws;
    
% Draw log-normal productivities of home and foreign firms
phiH = ones(max(MH),S)*epsilon;
phiF = ones(max(MF),S)*epsilon;
    
indexH = false(max(MH),S);
indexF = false(max(MF),S);
    
rng(aseed);
for j = 1:S
    phiH(1:MH(j),j) = exp(mu_H(j)+theta*randn(MH(j),1));
    indexH(1:MH(j),j) = true(MH(j),1);
    phiF(1:MF(j),j) = exp(theta*randn(MF(j),1));
    indexF(1:MF(j),j) = true(MF(j),1);
end

ZHS = phiH(1:max(sum(indexH(:,small))),small);
ZFS = phiF(1:max(sum(indexF(:,small))),small);
ZHL = phiH(:,~small);
ZFL = phiF(:,~small);

RTS = exp(mu_H(small));
RTL = exp(mu_H(~small));
RT = exp(mu_H);
%% Compute GE variables
   
% Guess home and foreign output Y
Y0=123;
YF0=2*Y0;

% Run loops to solve model (step 3 of estimation procedure)
[iter,Y,YF,LF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);

%% Compute 15 standard (target) moments

[StandardMom]=Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF);   

%% Compute additional moments

% Pareto-bound
paretonb = 75;

% Multiple of sectors (Total # of sectors = 117*4*scale)
scale = 1; %22;

% Run PE model
[ALPHA,PHIFVEC,LAMBDAFVEC,KVEC,KFVEC,KHH,KFH,MH,RT,AdditionalMom]=PEmoments(sigma,F,tau,theta,muT,sigmaT,w,wF,Y,YF,vMU,BER,paretonb,scale);

% Compute Verification Moments for Robustness 2

PKM = sum(KHH./MH'==1)/(S*scale);
EKM = mean(KHH./MH');
Corra = corrcoef(ALPHA',KHH./MH');
CorrT = corrcoef(RT,KHH./MH');

VerificationMom = [PKM,EKM,Corra(1,2),CorrT(1,2)];
%% Save moments

save('Results/Moments_seed1_grid2','StandardMom','AdditionalMom','VerificationMom')

%% Generate Graphs

% Compute required input

NP = 10;    % Number of bins
LPCT = prctile(LAMBDAFVEC,[1/NP:1/NP:1-1/NP]*100)'; % Percentiles
VPCT = zeros(NP,5);     % Matrix for plotting
VPCT(1,1) = mean(LAMBDAFVEC(LAMBDAFVEC<LPCT(1))>1.33*PHIFVEC(LAMBDAFVEC<LPCT(1)));  % White bar
VPCT(1,2) = mean(LAMBDAFVEC(LAMBDAFVEC<LPCT(1))>1.5*PHIFVEC(LAMBDAFVEC<LPCT(1)));   % Blue bar
VPCT(1,3) = mean(LAMBDAFVEC(LAMBDAFVEC<LPCT(1))>2*PHIFVEC(LAMBDAFVEC<LPCT(1)));     % Red bar
VPCT(1,4) = mean(LAMBDAFVEC(LAMBDAFVEC<LPCT(1))-PHIFVEC(LAMBDAFVEC<LPCT(1)));       % Blue bar (b)
VPCT(1,5) = sum(ALPHA(LAMBDAFVEC<LPCT(1))'.*LAMBDAFVEC(LAMBDAFVEC<LPCT(1)));         % Red bar (b)
% VPCT(1,:) = zeros(1,5);
for ip=2:NP-1
    VPCT(ip,1) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip))>1.33*PHIFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)));
    VPCT(ip,2) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip))>1.5*PHIFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)));
    VPCT(ip,3) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip))>2*PHIFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)));
    VPCT(ip,4) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)) - PHIFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)));
    VPCT(ip,5) = sum(ALPHA(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip))'.*LAMBDAFVEC(LAMBDAFVEC>=LPCT(ip-1)&LAMBDAFVEC<LPCT(ip)));
end
VPCT(NP,1) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(NP-1))>1.33*PHIFVEC(LAMBDAFVEC>=LPCT(NP-1)));
VPCT(NP,2) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(NP-1))>1.5*PHIFVEC(LAMBDAFVEC>=LPCT(NP-1)));
VPCT(NP,3) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(NP-1))>2*PHIFVEC(LAMBDAFVEC>=LPCT(NP-1)));
VPCT(NP,4) = mean(LAMBDAFVEC(LAMBDAFVEC>=LPCT(NP-1))-PHIFVEC(LAMBDAFVEC>=LPCT(NP-1)));
VPCT(NP,5) = sum(ALPHA(LAMBDAFVEC>=LPCT(NP-1))'.*LAMBDAFVEC(LAMBDAFVEC>=LPCT(NP-1)));

VPCT(:,5) = VPCT(:,5)/sum(ALPHA'.*LAMBDAFVEC);



% Figure 1(a)

figure(1)
bar(VPCT(:,1),'w','FaceAlpha',0)%,'BarWidth', 1)
hold on
bar(VPCT(:,2),'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33) %'r','FaceAlpha',0.65)%,'BarWidth', 1)
bar(VPCT(:,3),'FaceColor','[0.900, 0.20, 0.05]') %'b','FaceAlpha',0.5)%,'BarWidth', 1)
lgd = legend('$1/4$ Granular','$1/3$ Granular', '$1/2$ Granular');
set(lgd,'box','off','location','northwest','interpreter','latex','fontsize',20)
% xlim([0.5, 10.55])
% ylim([0, .3])
xlabel('Deciles of sectors, by export intensity $\Lambda_z^\ast$','interpreter','latex','fontsize',20)
set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
set(gca,'xticklabel',[0 round(LPCT'*100)/100 1.00])
% set(gca,'ytick',[0.05:.05:.25])
saveas(gcf,'Results/Graphs/Figure_1a','epsc')

% Figure 1(b)


figure(2)
bar(VPCT(:,5),'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33,'linestyle','--')
hold on
bar(VPCT(:,4),'FaceColor','[0.900, 0.20, 0.05]')
lgd = legend('Export share','Granular exports');
set(lgd,'box','off','location','northwest','interpreter','latex','fontsize',20)
xlim([0.5, 10.55])
xlabel('Deciles of sectors, by export intensity $\Lambda_z^\ast$','interpreter','latex','fontsize',20)
set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
set(gca,'xticklabel',[0 round(LPCT'*100)/100 1.00])
% set(gca,'ytick',[-0.05 0:.1:.3])
saveas(gcf,'Results/Graphs/Figure_1b','epsc')

% Generate the two verification histograms

figure(3)
hold on 
% edges = [0:20:900];
h = histogram(KHH./MH',20);
title('Distribution of $\tilde{K}_z/M_z$','interpreter','latex')
xlabel('$\tilde{K}_z/M_z$','interpreter','latex')
ylabel('Frequency')
saveas(gcf,'Results/Graphs/hist_KM.png')
hold off


figure(4)
hold on 
% edges = [0:20:900];
h = histogram(KFH./MF',20);
title('Distribution of $\tilde{K}_z^*/M_z^*$','interpreter','latex')
xlabel('$\tilde{K}_z^*/M_z^*$','interpreter','latex')
ylabel('Frequency')
saveas(gcf,'Results/Graphs/hist_KMstar.png')
hold off