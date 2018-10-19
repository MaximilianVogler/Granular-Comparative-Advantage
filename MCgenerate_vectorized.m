function [ ALPHA,RT,ZHS,ZFS,ZHL,ZFL ] = MCgenerate_vectorized(S,MH_small,MF_small,MH_large,MF_large,theta,muT,sigmaT)
% This code generates productivity parameters as in step 1 and 2 of the
% algorithm.
split_param=1.25;   % This determines the split between small and large sectors.
RT=exp(muT+sigmaT*randn(1,S));      % Make a draw for the relative fundamental comparative advantages from a lognormal distribution for each sector.

cdshares = csvread('cdshares.csv'); % Cobb-Douglas shares from external data source

ALPHA=cdshares;     % Rename CD shares in a vector
small=(ALPHA<split_param/S);    % Find indicator for small sectors (why do we have this exact cutoff rule?)
Nsmall=sum(small);  % Number of small sectors
Nlarge=S-Nsmall;    % Number of large sectors

RTS=RT(small);  % Relative productivity vector for small sectors
RTL=RT(~small);     % Relative productivity vector for large sectors

UH0S=exprnd(1,MH_small,Nsmall);  % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1
UF0S=exprnd(1,MF_small,Nsmall);  % Draw U of most productive small foreign shadow firm and spacings in each sector from exponential with mean 1
UHS=cumsum(UH0S);   % Cumulate to get U of each home firm in each sector
UFS=cumsum(UF0S);   % Cumulate to get U of each foreign firm in each sector
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1/theta);  % Transform into phi for home using T_z
ZFS=UFS.^(-1/theta);    % Transform into phi for foreign using T_Z^{ast}=1

UH0L=exprnd(1,MH_large,Nlarge); % Same as above for large shadow firms
UF0L=exprnd(1,MF_large,Nlarge);
UHL=cumsum(UH0L);
UFL=cumsum(UF0L);
ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1/theta);
ZFL=UFL.^(-1/theta);

end

