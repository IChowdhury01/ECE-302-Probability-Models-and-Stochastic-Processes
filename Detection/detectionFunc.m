function [Pf,Pm,Pd,Perror] = detectionFunc(A,X_mean,X_var,X_stddev,nTrials,p0,p1,sn,R)

X = normrnd(X_mean,X_stddev,1,nTrials);
Y = [X(1:p0*nTrials) (A + X((p0*nTrials) + 1:end))];
samp = sort(Y);

PH1 = normpdf(samp,A,X_stddev);
PH0 = normpdf(samp,X_mean,X_stddev);
n = PH1/PH0;

idx = min(find(samp > R));
Pf_e = 4*sum(PH0(idx:end))/sum(4*PH0);
Pm_e = sum(PH1(1:idx-1))/sum(PH1);
Pd_e = sum(PH1(idx:end))/sum(PH1);
    
Perror_e = p0*Pf_e + p1*Pm_e;

Pf = Pf_e;
Pm = Pm_e;
Pd = Pd_e;
Perror = Perror_e;

end