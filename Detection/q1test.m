%% ToDo
% Finish parts D and E
% 8.52) Plot theoretical curve alongside actual for figure 1
% Comment up code
% Create Publisher report
% Email submission

% Explain code (next week) 
%% Ivan Chowdhury, Josh Go
% ECE302 Probability Models & Stochastic Processes
% April 18, 2019
% Detection Project

clc;
clear all;

%% Question 1 - Radar Detection
% a)
% Givens
A = 1;
X_mean = 0;
X_var = 1;
X_stddev = sqrt(X_var);
nTrials = 1000;
p0 = 0.8;
p1 = (1-p0); %0.2

sn = A/(X_var);
R = A/2 + X_var*log(p0/p1)/A;

for i=1:nTrials
    [Pf,Pm,Pd,Perror_a] = detectionFunc(A,X_mean,X_var,X_stddev,nTrials,p0,p1,sn,R);
    Pf_e(:,i) = Pf;
    Pm_e(:,i) = Pm;
    Perror_E(:,i) = Perror_a;
end

for i = 1:nTrials
    detect = @(x) exp(-((x-A).^2 )/(2*X_var.^2))/((2*pi*X_var)^0.5);
    Pd = integral(detect,R,Inf);
    false = @(x) exp(-((x).^2 )/(2*X_var.^2))/((2*pi*X_var)^0.5);
    Pf = integral(false,R,Inf);
    Perror_t = p0*Pf + p1*(1-Pd);
    Perror_T(:,i) = Perror_t; 
end

PERRe = mean(Perror_E)
PERRt = mean(Perror_T)

%b)
inc = 1; %Increment
% A=1
A = 1;
sn = A/(X_var);
for i = 2:inc:1000
    [PfA1,PmA1,PdA1,~] = detectionFunc(A,X_mean,X_var,X_stddev,nTrials,i/1000,(1-(i/1000)),sn,(A/2 + X_var*log((i/1000)/(1-(i/1000)))/A));
    Pf_A1(i) = PfA1;
    Pd_A1(i) = PdA1;
end
% A=2
A = 2;
sn = A/(X_var); 
R = A/2 + X_var*log(p0/p1)/A;
for i = 2:inc:1000
    [PfA2,PmA2,PdA2,~] = detectionFunc(A,X_mean,X_var,X_stddev,nTrials,i/1000,(1-(i/1000)),sn,(A/2 + X_var*log((i/1000)/(1-(i/1000)))/A));
    Pf_A2(i) = PfA2;
    Pd_A2(i) = PdA2;
end
% A=3
A = 2.5;
sn = A/(X_var);
R = A/2 + X_var*log(p0/p1)/A;
for i = 2:inc:1000
    [PfA3,PmA3,PdA3,~] = detectionFunc(A,X_mean,X_var,X_stddev,nTrials,i/1000,(1-(i/1000)),sn,(A/2 + X_var*log((i/1000)/(1-(i/1000)))/A));
    Pf_A3(i) = PfA3;
    Pd_A3(i) = PdA3;
end


Pf_A1 = sort(Pf_A1); 
Pd_A1 = sort(Pd_A1);
Pf_A2 = sort(Pf_A2); 
Pd_A2 = sort(Pd_A2);
Pf_A3 = sort(Pf_A3); 
Pd_A3 = sort(Pd_A3);

figure;
hold on
plot(Pf_A1,Pd_A1)
xlabel('Pf')
ylabel('Pd')
plot(Pf_A2,Pd_A2)
xlabel('Pf')
ylabel('Pd')
plot(Pf_A3,Pd_A3)
title('Receiver Operating Curves for Different SNR')
xlabel('Pf')
ylabel('Pd')
legend('A=1','A=2','A=3')
hold off
%% c)
% Cost Structure
C01 = 10;   % Miss
C10 = 1;    % False Alarm
C00 = 0;    % Correct Reject 
C11 = 0;    % Detection

minmaxeq = Pd_A1 + Pf_A1./10;
minRiskPoint = min(find(minmaxeq > 0.1))    % Index used to find minimum risk decision rule

% Plot data
figure;
hold on
plot(Pf_A1,Pd_A1) %P(false alarm) vs P(detection), ranging from 0 to 1
plot(Pf_A1(minRiskPoint),Pd_A1(minRiskPoint),'*')
title('ROC Curve with minimum risk point marked')
xlabel('P(False Alarm)')
ylabel('P(Detection)')

%% d)

%% e)
