function [ K,KF,PHI,PHIF,LAMBDA,LAMBDAF,MU,MUF ] = PEreplication_vectorized( sigma,theta,F,tau,ALPHA,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,vMU )
% Takes parameters, random draws, and initial guess as inputs. Outputs are
% K, Phi, Lambda, and aggregate markup Mu in both countries. This is the
% matrix version of PEreplication.

MF_large=size(ZFL,1);
MF_small=size(ZFS,1);
MH_large=size(ZHL,1);
MH_small=size(ZHS,1); % get number of firms
S=size(ZHL,2)+size(ZHS,2); % fix number of sectors

split_param=1.25; % fix parameter that determines whether sectors are small
tol=1e-2; % set tolerance level for A-B loop

PHIHVEC=zeros(S,1);
MUHVEC=zeros(S,1);
PHIFVEC=zeros(S,1);
MUFVEC=zeros(S,1);
KVEC=zeros(S,1);
KFVEC=zeros(S,1);
LAMBDAHVEC=zeros(S,1);
LAMBDAFVEC=zeros(S,1); % initialize vectors that contain PE data

%% INITIALIZE DATA
small=(ALPHA<split_param/S); % identify small sectors

ALPHAS=ALPHA(small)'; 
ALPHAL=ALPHA(~small)';
RTS=RT(small);
RTL=RT(~small); % generate vectors of alphas/comparative advantages for large/small sectors

MCHS=w./ZHS;
MCFS=wF./ZFS; % small sector marginal costs

MCHL=w./ZHL;
MCFL=wF./ZFL; % large sector marginal costs

SCD=468;
FN=F*SCD/S; % adjust fixed cost of entry

%% SMALL SECTORS
MCHMS=[MCHS; tau*MCFS];
MCFMS=[MCFS; tau*MCHS]; 

[MCHMS,permH]=sort(MCHMS,1);
[MCFMS,permF]=sort(MCFMS,1); % sort marginal costs in home and foreign
IOTAH=(permH>MH_small);
IOTAF=(permF>MF_small); % find foreign/home entrants in home/foreign markets

MFMH_small=size(MCHMS,1); 

% set home markups/shares/entry in CMU case
MKP0H=sigma/(sigma-1); 
PHM=(MKP0H*MCHMS).^(1-sigma);
SLASTH=PHM./cumsum(PHM); % share of last entrant
checkmatH=(SLASTH>sigma/Y0*w*FN./(ones(MFMH_small,1)*ALPHAS)); % checks which firms enter in CMU case 
SHM=(PHM.*checkmatH)./(ones(MFMH_small,1)*sum(PHM.*checkmatH));
EPSH=sigma; % used to compute markups at the end of the loop

% same as above in foreign market
MKP0F=sigma/(sigma-1);
PFM=(MKP0F*MCFMS).^(1-sigma);
SLASTF=PFM./cumsum(PFM);
checkmatF=(SLASTF>sigma/YF0*wF*FN./(ones(MFMH_small,1)*ALPHAS));
SFM=(PFM.*checkmatF)./(ones(MFMH_small,1)*sum(PFM.*checkmatF));
EPSF=sigma;

% enter loop in VMU case
if vMU==1
    EPSH=sigma*(1-SHM)+SHM;
    EPSH(EPSH<=1.25)=1.25; % cap markup EPSH/(EPSH-1) at 5
    
    EPSF=sigma*(1-SFM)+SFM;
    EPSF(EPSF<=1.25)=1.25; % same in foreign
    
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;
        
        % update shares in home market
        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1)); % set markups for current step
        PHM=(MKPNH.*MCHMS).^(1-sigma); % set prices given markups
        SHM=(PHM.*checkmatH)./(ones(MFMH_small,1)*sum(PHM.*checkmatH)); % calculate shares given prices
        EPSH=sigma*(1-SHM)+SHM; % update next step epsilon
        EPSH(EPSH<=1.25)=1.25; % cap markups at 5 for next step
        
        % same in foreign market       
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFMS).^(1-sigma);
        SFM=(PFM.*checkmatF)./(ones(MFMH_small,1)*sum(PFM.*checkmatF));
        EPSF=sigma*(1-SFM)+SFM;
        EPSF(EPSF<=1.25)=1.25;
        
        % check for convergence
        diff=sum(sum(abs([MKPNH-MKP0H;MKPNF-MKP0F])));
        
        % update markups
        MKP0H=MKPNH;
        MKP0F=MKPNF;

    end

    % check entry conditions again with variable markups
    PHM=(MKP0H.*MCHMS).^(1-sigma);
    SLASTH=PHM./cumsum(PHM); % share of last entrant
    checkmatH=(SLASTH>sigma/Y0*w*FN./(ones(MFMH_small,1)*ALPHAS)); % determine which firms enter with CMU fringe
    SHM=(PHM.*checkmatH)./(ones(MFMH_small,1)*sum(PHM.*checkmatH)); % update shares
    EPSH=sigma*(1-SHM)+SHM;
    EPSH(EPSH<=1.25)=1.25; % cap markups at 5 (for later)
    
    % same in foreign market
    PFM=(MKP0F.*MCFMS).^(1-sigma);
    SLASTF=PFM./cumsum(PFM);
    checkmatF=(SLASTF>sigma/YF0*wF*FN./(ones(MFMH_small,1)*ALPHAS));
    SFM=(PFM.*checkmatF)./(ones(MFMH_small,1)*sum(PFM.*checkmatF));
    EPSF=sigma*(1-SFM)+SFM;
    EPSF(EPSF<=1.25)=1.25;
    
    % same loop as before, but with new entrants
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHMS).^(1-sigma);
        SHM=(PHM.*checkmatH)./(ones(MFMH_small,1)*sum(PHM.*checkmatH));

        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFMS).^(1-sigma);
        SFM=(PFM.*checkmatF)./(ones(MFMH_small,1)*sum(PFM.*checkmatF));

        EPSH=sigma*(1-SHM)+SHM;
        EPSH(EPSH<=1.25)=1.25;

        EPSF=sigma*(1-SFM)+SFM;
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

KVEC(small)=sum(checkmatH); % vector of number of firms that enter in home
KFVEC(small)=sum(checkmatF);
PHIHVEC(small)=1./(1+(tau*wF/w)^theta.*RTS); % expected import share vector in home
PHIFVEC(small)=1./(1+(tau*w/wF)^theta./RTS);
MUHVEC(small)=sum(SHM./MUHM); % average inverse markup vector in home
MUFVEC(small)=sum(SFM./MUFM);
LAMBDAHVEC(small)=sum(IOTAH.*SHM); % realized import share vector in home
LAMBDAFVEC(small)=sum(IOTAF.*SFM);

%% LARGE SECTORS
% SAME CODE AS ABOVE
MCHML=[MCHL; tau*MCFL];
MCFML=[MCFL; tau*MCHL];

[MCHML,permH]=sort(MCHML,1);
[MCFML,permF]=sort(MCFML,1);
IOTAH=(permH>MH_large);
IOTAF=(permF>MF_large);

MFMH_large=size(MCHML,1);

MKP0H=sigma/(sigma-1);
PHM=(MKP0H*MCHML).^(1-sigma);
SLASTH=PHM./cumsum(PHM);
checkmatH=(SLASTH>sigma/Y0*w*FN./(ones(MFMH_large,1)*ALPHAL));
SHM=(PHM.*checkmatH)./(ones(MFMH_large,1)*sum(PHM.*checkmatH));
EPSH=sigma;

MKP0F=sigma/(sigma-1);
PFM=(MKP0F*MCFML).^(1-sigma);
SLASTF=PFM./cumsum(PFM);
checkmatF=(SLASTF>sigma/YF0*wF*FN./(ones(MFMH_large,1)*ALPHAL));
SFM=(PFM.*checkmatF)./(ones(MFMH_large,1)*sum(PFM.*checkmatF));
EPSF=sigma;

if vMU==1
    EPSH=sigma*(1-SHM)+SHM;
    EPSH(EPSH<=1.25)=1.25;
    
    EPSF=sigma*(1-SFM)+SFM;
    EPSF(EPSF<=1.25)=1.25;
    
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHML).^(1-sigma);
        SHM=(PHM.*checkmatH)./(ones(MFMH_large,1)*sum(PHM.*checkmatH));

        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFML).^(1-sigma);
        SFM=(PFM.*checkmatF)./(ones(MFMH_large,1)*sum(PFM.*checkmatF));

        EPSH=sigma*(1-SHM)+SHM;
        EPSH(EPSH<=1.25)=1.25;

        EPSF=sigma*(1-SFM)+SFM;
        EPSF(EPSF<=1.25)=1.25;

        diff=sum(sum(abs([MKPNH-MKP0H;MKPNF-MKP0F])));

        MKP0H=MKPNH;
        MKP0F=MKPNF;

    end

    PHM=(MKP0H.*MCHML).^(1-sigma);
    SLASTH=PHM./cumsum(PHM);
    checkmatH=(SLASTH>sigma/Y0*w*FN./(ones(MFMH_large,1)*ALPHAL));

    SHM=(PHM.*checkmatH)./(ones(MFMH_large,1)*sum(PHM.*checkmatH));
    EPSH=sigma*(1-SHM)+SHM;
    EPSH(EPSH<=1.25)=1.25;

    PFM=(MKP0F.*MCFML).^(1-sigma);
    SLASTF=PFM./cumsum(PFM);
    checkmatF=(SLASTF>sigma/YF0*wF*FN./(ones(MFMH_large,1)*ALPHAL));

    SFM=(PFM.*checkmatF)./(ones(MFMH_large,1)*sum(PFM.*checkmatF));
    EPSF=sigma*(1-SFM)+SFM;
    EPSF(EPSF<=1.25)=1.25;

    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHML).^(1-sigma);
        SHM=(PHM.*checkmatH)./(ones(MFMH_large,1)*sum(PHM.*checkmatH));

        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFML).^(1-sigma);
        SFM=(PFM.*checkmatF)./(ones(MFMH_large,1)*sum(PFM.*checkmatF));

        EPSH=sigma*(1-SHM)+SHM;
        EPSH(EPSH<=1.25)=1.25;

        EPSF=sigma*(1-SFM)+SFM;
        EPSF(EPSF<=1.25)=1.25;

        diff=sum(sum(abs([MKPNH-MKP0H;MKPNF-MKP0F])));

        MKP0H=MKPNH;
        MKP0F=MKPNF;

    end
end

checkmatH=(SHM>0);
checkmatF=(SFM>0);
MUHM=EPSH./(EPSH-1);
MUFM=EPSF./(EPSF-1);

KVEC(~small)=sum(checkmatH);
KFVEC(~small)=sum(checkmatF);
PHIHVEC(~small)=1./(1+(tau*wF/w)^theta.*RTL);
PHIFVEC(~small)=1./(1+(tau*w/wF)^theta./RTL);
MUHVEC(~small)=sum(SHM./MUHM);
MUFVEC(~small)=sum(SFM./MUFM);
LAMBDAHVEC(~small)=sum(IOTAH.*SHM);
LAMBDAFVEC(~small)=sum(IOTAF.*SFM);

% Take averages of vectors to obtain aggregates
K=mean(KVEC); % number of firms in home
KF=mean(KFVEC);
PHI=mean(PHIHVEC); % expected import share in home
PHIF=mean(PHIFVEC);
LAMBDA=mean(LAMBDAHVEC); % realized import share in home
LAMBDAF=mean(LAMBDAFVEC);
MU=1./sum(ALPHA.*MUHVEC); % aggregate markup in home
MUF=1./sum(ALPHA.*MUFVEC);

end


