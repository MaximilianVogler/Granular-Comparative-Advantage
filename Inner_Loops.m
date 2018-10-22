function [KVEC,KFVEC,PHIHVEC,PHIFVEC,MUHVEC,MUFVEC,LAMBDAHVEC,LAMBDAFVEC] = Inner_Loops(sigma,theta,F,tau,ALPHA,RT,ZH,ZF,w,wF,Y0,YF0,vMU,BER)

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
PHM=(MKP0H*MCHM).^(1-sigma);                                                % Why is there no price index???
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
        PHM=(MKPNH.*MCHM).^(1-sigma); % set prices given markups            % PRICE INDEX???
        % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH)); % calculate shares given prices
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         % TEST!!!
        EPSH=vMU_markup(sigma,SHM,BER);
        EPSH(EPSH<=1.25)=1.25; % cap markups at 5 for next step
        
        % same in foreign market       
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                       % PRICE INDEX???
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
    PHM=(MKP0H.*MCHM).^(1-sigma);                                           % PRICE INDEX???
    SLASTH=PHM./cumsum(PHM); % share of last entrant
    % checkmatH=(SLASTH>sigma/Y0*w*F./(ones(MFMH,1)*ALPHA)); % determine which firms enter with CMU fringe
    checkmatH=(SLASTH>sigma/Y0*w*F./(repmat(ALPHA,MFMH,1)));                % TEST!!!
    % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH)); % update shares
    SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));             % TEST!!!
    EPSH=vMU_markup(sigma,SHM,BER);
    EPSH(EPSH<=1.25)=1.25; % cap markups at 5 (for later)
    
    % same in foreign market
    PFM=(MKP0F.*MCFM).^(1-sigma);                                           % PRICE INDEX???
    SLASTF=PFM./cumsum(PFM);
    % checkmatF=(SLASTF>sigma/YF0*wF*F./(ones(MFMH,1)*ALPHA));      
    checkmatF=(SLASTF>sigma/YF0*wF*F./(repmat(ALPHA,MFMH,1)));              % TEST!!!
    % SFM=(PFM.*checkmatF)./(ones(MFMH,1)*sum(PFM.*checkmatF));
    SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMs,1));             % TEST!!!
    EPSF=vMU_markup(sigma,SFM,BER);
    EPSF(EPSF<=1.25)=1.25;
    
    % same loop as before, but with new entrants
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHM).^(1-sigma);                                       % PRICE INDEX???
        % SHM=(PHM.*checkmatH)./(ones(MFMH,1)*sum(PHM.*checkmatH));
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         % TEST!!!
        
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                       % PRICE INDEX???  
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
MUHVEC=sum(SHM./MUHM); % average inverse markup vector in home
MUFVEC=sum(SFM./MUFM);
LAMBDAHVEC=sum(IOTAH.*SHM); % realized import share vector in home
LAMBDAFVEC=sum(IOTAF.*SFM);

end