function [mom] = Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF)
% This function computes the 15 target moments required for estimation.

mom(1:2) = moment_stats(log(KHH(KHH>0)));       % Mean and standard deviation of home firms active in home.

%DSHM = sort(DSHM,'descend');
%DSHM = DSHM./repmat(sum(DSHM),size(DSHM,1),1); % Divide to make the share relative     NEED TO TEST!!!
%TOP1 = DSHM(1,:);
mom(3:4) = moment_stats(TOP1(TOP1>0));

%TOP3 = sum(DSHM(1:3,:));
mom(5:6) = moment_stats(TOP3(TOP3>0));

mom(7:8) = moment_stats(LAMBDAHVEC(LAMBDAHVEC>0));

% Total Exports for each sector
X = ALPHA'.*YF.*XS;

% Domestic sales for each sector
YX = ALPHA'.*Y.*YXS;

LAMBDAPRIME = LAMBDAFVEC*YF/Y;
mom(9:10) = moment_stats(LAMBDAPRIME(LAMBDAPRIME>0));
% Verification:
%XXX = X./YX./YXS;
%mom(9:10) = moment_stats(XXX(XXX>0));

mom(11) = mean(X>YX);

% Regress Lambda on TOP1

Control = log(1+ALPHA'.*Y.*(1-LAMBDAHVEC));

stats = regstats(LAMBDAPRIME',[TOP1', Control'],'linear',{'beta','covb'});
mom(12) = stats.beta(2);

stats = regstats(LAMBDAPRIME',[TOP3', Control'],'linear',{'beta','covb'});
mom(13) = stats.beta(2);

stats = regstats(LAMBDAHVEC,[TOP1', Control'],'linear',{'beta','covb'});
mom(14) = stats.beta(2);

stats = regstats(LAMBDAHVEC,[TOP3', Control'],'linear',{'beta','covb'});
mom(15) = stats.beta(2);

end