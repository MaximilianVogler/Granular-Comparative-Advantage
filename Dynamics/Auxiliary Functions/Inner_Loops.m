function [KVEC,KFVEC,PHIHVEC,PHIFVEC,MUHVEC,MUFVEC,LAMBDAHVEC,LAMBDAFVEC,KHH,TOP1,TOP3,XS,YXS,Paretovec,varphi_bar,DVEC,DSHM1,DSHM] = Inner_Loops(sigma,theta,F,tau,ALPHA,RT,ZH,ZF,w,wF,Y0,YF0,vMU,BER,paretonb,AddMom,S)

tol=1e-2; % set tolerance level for A-B loop

%% Set up loops
MH=size(ZH,1);
MF=size(ZF,1);

MCH=w./ZH;
MCF=wF./ZF;

F = F*476/S;                                                              % This is required to adjust the fixed costs for the simulation size. This is because we have chosen the normalization sum(alpha)=1.  
MCHM=[MCH; tau*MCF];                                                      % This implies that the size of the sector changes as the average alpha decreases as the number of sectors S increases. 
MCFM=[MCF; tau*MCH];                                                      % To keep the model stable, we thus need to adjust the fixed cost F accordingly. Equivalently, we could have had sum(alpha)=S, since this would have implied a constant average alpha of 1.

[MCHM,permH]=sort(MCHM,1);
[MCFM,permF]=sort(MCFM,1); % sort marginal costs in home and foreign
IOTAH=(permH>MH);
IOTAF=(permF>MF); % find foreign/home entrants in home/foreign markets

MFMH=size(MCHM,1); 

% JJ are the indices that transform the sorted matrix back into the
% original matrix
% [~,JJ] = sort(permH);

%% Constant Markup
% Set home markups/shares/entry in CMU case
MKP0H=sigma/(sigma-1); 
PHM=(MKP0H*MCHM).^(1-sigma);                                                
SLASTH=PHM./cumsum(PHM); % share of last entrant
checkmatH=(SLASTH>sigma/Y0*w*F./(repmat(ALPHA,MFMH,1)));                   % checks which firms enter in CMU case 
SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));                
EPSH=sigma; % used to compute markups at the end of the loop

% same as above in foreign market
MKP0F=sigma/(sigma-1);
PFM=(MKP0F*MCFM).^(1-sigma);
SLASTF=PFM./cumsum(PFM);     
checkmatF=(SLASTF>sigma/YF0*wF*F./(repmat(ALPHA,MFMH,1)));                 
SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));                  
EPSF=sigma;

%% Variable Markup 
% enter loop in VMU case
if vMU==1
    EPSH=vMU_markup(sigma,SHM,BER);
    EPSH(EPSH<=1.25)=1.25; % cap markup EPSH/(EPSH-1) at 5                  
    
    EPSF=vMU_markup(sigma,SFM,BER);
    EPSF(EPSF<=1.25)=1.25; % same in foreign
    
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;
        
        % update shares in home market
        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1)); % set markups for current step
        PHM=(MKPNH.*MCHM).^(1-sigma); % set prices given markups            
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         % calculate shares given prices
        EPSH=vMU_markup(sigma,SHM,BER);
        EPSH(EPSH<=1.25)=1.25; % cap markups at 5 for next step
        
        % same in foreign market       
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                       
        SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));         
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
    checkmatH=(SLASTH>sigma/Y0*w*F./(repmat(ALPHA,MFMH,1)));                % determine which firms enter with CMU fringe  
    SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));             % update shares
    EPSH=vMU_markup(sigma,SHM,BER);
    EPSH(EPSH<=1.25)=1.25; % cap markups at 5 (for later)
    
    % same in foreign market
    PFM=(MKP0F.*MCFM).^(1-sigma);                                          
    SLASTF=PFM./cumsum(PFM);     
    checkmatF=(SLASTF>sigma/YF0*wF*F./(repmat(ALPHA,MFMH,1)));              
    SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));          
    EPSF=vMU_markup(sigma,SFM,BER);
    EPSF(EPSF<=1.25)=1.25;
    
    % same loop as before, but with new entrants
    iter=0;
    diff=1;
    while (diff>tol) && (iter<101)
        iter=iter+1;

        MKPNH=1/2*(MKP0H+EPSH./(EPSH-1));
        PHM=(MKPNH.*MCHM).^(1-sigma);                                     
        SHM=(PHM.*checkmatH)./(repmat(sum(PHM.*checkmatH),MFMH,1));         
        
        MKPNF=1/2*(MKP0F+EPSF./(EPSF-1));
        PFM=(MKPNF.*MCFM).^(1-sigma);                                           
        SFM=(PFM.*checkmatF)./(repmat(sum(PFM.*checkmatF),MFMH,1));         
        
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
% MUHVEC=sum(SHM./MUHM); % average inverse markup vector in home
% MUFVEC=sum(SFM./MUFM);
MUHVEC=sum((1-IOTAH).*SHM./MUHM);
MUFVEC=sum((1-IOTAF).*SFM./MUFM);                                           
LAMBDAHVEC=sum(IOTAH.*SHM); % realized import share vector in home
LAMBDAFVEC=sum(IOTAF.*SFM);
DVEC=sum((1-IOTAH).*SHM);   % domestic sales of domestic firms

% Auxiliary variables for moments
KHH = sum(checkmatH.*(1-IOTAH));    % Number of home firms active in home for each sector
KHF = sum(checkmatF.*IOTAF);        % Number of home firms active in foreign for each sector
KFH = sum(checkmatH.*IOTAH);        % Number of foreign firms active in home for each sector
KFF = sum(checkmatF.*(1-IOTAF));    % Number of foreign firms active in foreign for each sector
DSHM = SHM.*checkmatH.*(1-IOTAH);   % Share on the home market relative to other domestic firms (equation 17)
% DSHM1 = DSHM(JJ);
DSHM1 = sort_back(DSHM,permH);
DSHM1 = DSHM1./repmat(sum(DSHM1),size(DSHM1,1),1);
DSHM = sort(DSHM,'descend');
DSHM = DSHM./repmat(sum(DSHM),size(DSHM,1),1); % Divide to make the share relative    
TOP1 = DSHM(1,:);
TOP3 = sum(DSHM(1:3,:));
XS = sum(IOTAF.*SFM);
YXS = 1-sum(IOTAH.*SHM);

N = length(DSHM(1,:));
Paretovec=zeros(1,N);

if AddMom == 1
    Y=log((1:MFMH)-0.5)';
    for i=1:N
    lshare=log(DSHM(:,i)); % generate log shares
    X=[0*Y+1 lshare];
    pareto_thresh=prctile(lshare(isfinite(lshare)),paretonb);
    index=isfinite(lshare)&lshare>pareto_thresh; % include if share is nonzero and sufficiently large
    if sum(index)>10 % exclude sectors with few firms
        B=regress(Y(index),X(index,:));
        Paretovec(i)=-B(2);
    else
        Paretovec(i)=0; 
    end
    end
end

varphi_bar = [KHH;KHF;KFH;KFF];
end