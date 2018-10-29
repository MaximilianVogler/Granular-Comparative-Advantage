% Testing Inner Loops

clear all;
close all;
clc;

sigma = 5;
theta = 4.306666666666667;
F = 0.945561052525253*10^(-5);
tau = 1.340707070707071;


aseed=1;

MF_large=10000;
MF_small=1400;
MH_large=5000;
MH_small=700;
% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4;                                         
S_init = length(cdshares_init);                         % Number of sectors in data
S = S_init*S_multiple;                                  % Number of sectors used (see footnote 56)

% Assign CD-shares across sectors 
ALPHA = cdshares_init;      

for iloop = 1:S_multiple-1;
    ALPHA = [ALPHA;cdshares_init(randperm(S_init))];
end
ALPHA = ALPHA/S_multiple;

% Split sectors into large and small based on CD-share.
split_param = 1.25;
small = (ALPHA<split_param/S);                            
Nsmall = sum(small);                           
Nlarge = S-Nsmall;

ALPHA = ALPHA(small)';

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

% Fix sigma
SIGMA = 5;

% Labor Normalization
L0 = 100;

% Input parameters
muT = 0.136767676767677;
sigmaT = 1.421717171717172;

RT=exp(muT+sigmaT*rtdraws);
RT=RT(small);
RTS=RT;

    ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
    ZFS=UFS.^(-1/theta);
    
    ZH=ZHS;
    ZF=ZFS;
    
    
    wR = 1.13;      
w = 1;          
wF = 1/wR; 


BER = 1;

vMU=1;
Y0 = 126.7506707854941;
YF0 = 1.526*Y0;

tol=1e-2; % set tolerance level for A-B loop

%% Set up loops
MH=size(ZH,1);
MF=size(ZF,1);

MCH=w./ZH;
MCF=wF./ZF;

MCHM=[MCH; tau*MCF];
MCFM=[MCF; tau*MCH]; 

[MCHM,permH]=sort(MCHM,1);
[MCFM,permF]=sort(MCFM,1); % sort marginal costs in home and foreign
IOTAH=(permH>MH);
IOTAF=(permF>MF); % find foreign/home entrants in home/foreign markets

MFMH=size(MCHM,1); 

%% Constant Markup
% Set home markups/shares/entry in CMU case
MKP0H=sigma/(sigma-1); 
PHM=(MKP0H*MCHM).^(1-sigma);                                                
SLASTH=PHM./cumsum(PHM); % share of last entrant
% checkmatH=(SLASTH>sigma/Y0*w*F./(ones(MFMH,1)*ALPHA)); % checks which firms enter in CMU case 
checkmatH=(SLASTH>sigma/Y0*w*F./(repmat(ALPHA,MFMH,1)));                    % TEST!!!
% SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH));
SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));                 % TEST!!!
EPSH=sigma; % used to compute markups at the end of the loop

% same as above in foreign market
MKP0F=sigma/(sigma-1);
PFM=(MKP0F*MCFM).^(1-sigma);
SLASTF=PFM./cumsum(PFM);
% checkmatF=(SLASTF>sigma/YF0*wF*F./(ones(MFMH,1)*ALPHA));      
checkmatF=(SLASTF>sigma/YF0*wF*F./(repmat(ALPHA,MFMH,1)));                  % TEST!!!
% SFM=(PFM.*checkmatF)./(ones(MFMH,1)*sum(PFM.*checkmatF));
SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));                 % TEST!!!    
EPSF=sigma;

%% Variable Markup 
% enter loop in VMU case
if vMU==1
    EPSH=vMU_markup(sigma,SHM,BER);
    EPSH(EPSH<=1.25)=1.25; % cap markup EPSH/(EPSH-1) at 5                  % WHY???
    
    EPSF=vMU_markup(sigma,SFM,BER);
    EPSF(EPSF<=1.25)=1.25; % same in foreign
    
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;
        
        % update shares in home market
        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1)); % set markups for current step
        PHM=(MKPNH.*MCHM).^(1-sigma); % set prices given markups            
        % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH)); % calculate shares given prices
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         % TEST!!!
        EPSH=vMU_markup(sigma,SHM,BER);
        EPSH(EPSH<=1.25)=1.25; % cap markups at 5 for next step
        
        % same in foreign market       
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                       
        % SFM=(PFM.*checkmatF)./(ones(MFMH,1)*sum(PFM.*checkmatF));
        SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));         % TEST!!!
        EPSF=vMU_markup(sigma,SFM,BER);
        EPSF(EPSF<=1.25)=1.25;
        
        % check for convergence
        diff=sum(sum(abs([MKPNH-MKP0H;MKPNF-MKP0F])));
        
        % update markups
        MKP0H=MKPNH;
        MKP0F=MKPNF;

    end
    
    % check entry conditions again with variable markups
    PHM=(MKP0H.*MCHM).^(1-sigma);                                           
    SLASTH=PHM./cumsum(PHM); % share of last entrant
    % checkmatH=(SLASTH>sigma/Y0*w*F./(ones(MFMH,1)*ALPHA)); % determine which firms enter with CMU fringe
    checkmatH=(SLASTH>sigma/Y0*w*F./(repmat(ALPHA,MFMH,1)));                % TEST!!!
    % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH)); % update shares
    SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));             % TEST!!!
    EPSH=vMU_markup(sigma,SHM,BER);
    EPSH(EPSH<=1.25)=1.25; % cap markups at 5 (for later)
    
    % same in foreign market
    PFM=(MKP0F.*MCFM).^(1-sigma);                                          
    SLASTF=PFM./cumsum(PFM);
    % checkmatF=(SLASTF>sigma/YF0*wF*F./(ones(MFMH,1)*ALPHA));      
    checkmatF=(SLASTF>sigma/YF0*wF*F./(repmat(ALPHA,MFMH,1)));              % TEST!!!
    % SFM=(PFM.*checkmatF)./(ones(MFMH,1)*sum(PFM.*checkmatF));
    SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));             % TEST!!!
    EPSF=vMU_markup(sigma,SFM,BER);
    EPSF(EPSF<=1.25)=1.25;
    
    % same loop as before, but with new entrants
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHM).^(1-sigma);                                     
        % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH));
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         % TEST!!!
        
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                      
        % SFM=(PFM.*checkmatF)./(ones(MFMH,1)*sum(PFM.*checkmatF));     
        SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));         % TEST!!!
        
        EPSH=vMU_markup(sigma,SHM,BER);
        EPSH(EPSH<=1.25)=1.25;
        
        EPSF=vMU_markup(sigma,SFM,BER);
        EPSF(EPSF<=1.25)=1.25;

        diff=sum(sum(abs([MKPNH-MKP0H;MKPNF-MKP0F])));

        MKP0H=MKPNH;
        MKP0F=MKPNF;

    end
end

checkmatH=(SHM>0); % matrix of firms that enter in home in equilibrium above
checkmatF=(SFM>0); 
MUHM=EPSH./(EPSH-1); % markups of firms in home
MUFM=EPSF./(EPSF-1);

%% Compute Summary Statistics 
KVEC=sum(checkmatH); % vector of number of firms that enter in home
KFVEC=sum(checkmatF);
PHIHVEC=1./(1+(tau*wF/w)^theta.*RT); % expected import share vector in home
PHIFVEC=1./(1+(tau*w/wF)^theta./RT);
%MUHVEC=sum(SHM./MUHM); % average inverse markup vector in home
%MUFVEC=sum(SFM./MUFM);
MUHVEC=sum((1-IOTAH).*SHM./MUHM);
MUFVEC=sum((1-IOTAF).*SFM./MUFM);                                           % Careful!!! Need to check whether this is true.
LAMBDAHVEC=sum(IOTAH.*SHM); % realized import share vector in home
LAMBDAFVEC=sum(IOTAF.*SFM);

% Auxiliary variables for moments
KHH = sum(checkmatH.*(1-IOTAH));    % Number of home firms active in home for each sector
DSHM = SHM.*checkmatH.*(1-IOTAH);    % Share on the home market relative to other domestic firms (equation 17)
DSHM = sort(DSHM,'descend');
DSHM = DSHM./repmat(sum(DSHM),size(DSHM,1),1); % Divide to make the share relative     NEED TO TEST!!!
TOP1 = DSHM(1,:);
%TOP1 = TOP1(TOP1>0);
XS = sum(IOTAF.*SFM);
YXS = 1-sum(IOTAH.*SHM);
