function [mom] = Moments(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,L0,small,vMU,BER)
% This function computes the 15 target moments required for estimation.

KHH = sum(checkmatH.*(1-IOTAH));    % Number of home firms active in home for each sector
mom(1:2) = moment_stats(log(KHH));       % Mean and standard deviation of home firms active in home.

DSHM = SHM.*checkmatH.*(1-IOTAH);    % Share on the home market relative to other domestic firms (equation 17)
DSHM = sort(DSHM,'descend');
DSHM = DSHM./repmat(sum(DSHM),size(DSHM,1),1); % Divide to make the share really realtive NEED TO TEST!!!
TOP1 = DSHM(1,:);
mom(3:4) = moment_stats(TOP1);

TOP3 = sum(DSHM(1:3,:));
mom(5:6) = moment_stats(TOP3);

mom(7:8) = moment_stats(LAMBDAHVEC);

% Import share in foreign, i.e. share of home exports in foreign for each
% sector (vector)
XS = sum(IOTAF.*SFM);

% Total Exports for each sector
X = ALPHA'.*YF.*XS;

% Domestic sales for each sector
YX = ALPHA'.*Y.*(1-sum(IOTAH.*SHM));

LAMBDAPRIME = LAMBDAFVEC*Y/YF;
mom(9:10) = moment_stats(LAMBDAPRIME);
% Verification:
% mom(9:10) = moment_stats(X./YX);

mom(11) = mean(X>YX);

% Regress Lambda on TOP1

B = regress(LAMBDAHVEC',[TOP1' log(ALPHA*Y)]);

mom(12) = B(1);

B = regress(LAMBDAHVEC',[TOP3' log(ALPHA*Y)]);

mom(13) = B(1);

B = regress(LAMBDAPRIME', [TOP1' log(ALPHA*Y)]);

mom(14) = B(1);

B = regress(LAMBDAPRIMT', [TOP3' log(ALPHA*Y)]);

mom(15) = B(1);

end