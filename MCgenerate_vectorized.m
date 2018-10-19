function [ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL] = MCgenerate_vectorized(S_multiple,split_param,outer,MH_small,MF_small,MH_large,MF_large,theta,muT,sigmaT)
%function [ALPHA,RT] = MCgenerate_vectorized(S_multiple)
% This code generates productivity parameters as in step 1 and 2 of the
% algorithm.

if outer == 1                                           % This will be called outside the estimation loop
    
% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_init = length(cdshares_init);                         % Number of sectors in data
S = S_init*S_multiple;                                  % Number of sectors used (see footnote 56)

% Assign CD-shares across sectors
ALPHA = cdshares_init;      

for iloop = 1:S_multiple-1;
    ALPHA = [ALPHA;cdshares_init(randperm(S_init))];
end
ALPHA = ALPHA/S_multiple;

% Split sectors into large and small based on number of CD-share.
small=(ALPHA<split_param/S);                            
Nsmall=sum(small);                           
Nlarge=S-Nsmall;

end


    
if outer == 0                                           % This will be called inside the estimation loop
    
% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z
% (step 1 of estimation procedure)
RT = exp(muT+sigmaT*randn(1,S));  
RTS = RT(small);
RTL = RT(~small);


% 
% UH0S=exprnd(1,MH_small,Nsmall);                % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1
% UF0S=exprnd(1,MF_small,Nsmall);                % Draw U of most productive small foreign shadow firm and spacings in each sector from exponential with mean 1
% UHS=cumsum(UH0S);                              % Cumulate to get U of each home firm in each sector
% UFS=cumsum(UF0S);                              % Cumulate to get U of each foreign firm in each sector
% ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1/theta); % Transform into phi for home using T_z
% ZFS=UFS.^(-1/theta);                           % Transform into phi for foreign using T_Z^{ast}=1
% 
% UH0L=exprnd(1,MH_large,Nlarge);                % Same as above for large shadow firms
% UF0L=exprnd(1,MF_large,Nlarge);
% UHL=cumsum(UH0L);
% UFL=cumsum(UF0L);
% ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1/theta);
% ZFL=UFL.^(-1/theta);

end

