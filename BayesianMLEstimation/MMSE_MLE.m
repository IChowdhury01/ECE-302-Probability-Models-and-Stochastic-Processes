% Ivan Chowdhury, Josh Go
% ECE302: Probability and Stochastic Processes
% MATLAB Project 1: Bayesian MMSE and MLE Estimators
% 4/1/2019

clear all;
close all; 
clc;
%% Scenario 1

% 1(a) - Bayesian MMSE Estimator

% Simulation Parameters
numIter = 500;     % Number of iterations used for simulation
numMeas =  50;     % Max number of measurements taken

% Estimation Parameters
h = 0.5;            % Known parameter
mean_theta = 4;     % Mean of Theta
var_theta = 1;      % Variance of Theta
mean_v = 0;         % Mean of the gaussian random variable v
var_v = 4;          % Variance of v

% Pre-allocate memory
MSE_measurement = zeros(1,numMeas);         % Array containing mean square errors for all measurements
MSE_iteration = zeros(numIter,numMeas);     % Array containing mean square errors for all iterations
MSE_total = zeros(3,numMeas);               % Array containing total MSE for each tested distribution

for j = 1:3 % Chooses which distribution theta is drawn from (for part d)
    for n = 1:numIter
        x = zeros(1,numMeas);   % Create array to hold each measurement of x = h*theta + v
        theta_dist = [normrnd(mean_theta,sqrt(var_theta)) 2*mean_theta*rand exprnd(1/mean_theta)]; % Generates sample from Normal, Uniform, and Exponential distributions using theta
        theta = theta_dist(j); % Value of theta pulled from random variable distribution.
        
        for k = 1:numMeas
            x(k)= h*theta + normrnd(mean_v,sqrt(var_v)); % Compute x measurement
            x_bar = sum(x)/k; % Average of x measurements for this iteration
            theta_est = (k/var_v*x_bar+mean_theta/var_theta)/((k/var_v + 1/var_theta)*h); % Bayesian MMSE of theta calculation
            
            MSE_measurement(k) = (theta - theta_est)^2; % Mean square error of Baye's estimation for each measurement
        end
        MSE_iteration(n,:) = MSE_measurement;    % MSE of Baye's estimation for each iteration
    end
MSE_total(j,:) = mean(MSE_iteration); % MSE for all tested distributions
end

BayesMSE = MSE_total(1,:);              % Mean square errors with gaussian prior
BayesMSE_uniform = MSE_total(2,:);      % MSE with theta pulled from uniform distribution 
BayesMSE_exp = MSE_total(3,:);          % MSE with theta pulled from exponential distribution

% 1(b) - ML Estimator

% Pre-allocate memory
MSE_measurement2 = zeros(1,numMeas);   % Arrays used to store MSE data for measurements and iterations
MSE_iteration2 = zeros(numIter,numMeas);

for n = 1:numIter
    x = zeros(1,numMeas);   % Create array to store x = h*theta + v values
    for k = 1:numMeas
        x(k) = h*mean_theta + normrnd(mean_v,sqrt(var_v)); % x = h*theta + v, pulled from gaussian distribution
        x_bar = sum(x)/k; % mean of x measurements
        
        theta_est2 = x_bar/h; % ML estimator of theta calculation
        MSE_measurement2(k)= ((mean_theta)- theta_est2)^2 ; % Mean square error of our ML estimator
    end
    MSE_iteration2(n,:) = MSE_measurement2; % Store mean square error across all iterations
end

ML_MSE = mean(MSE_iteration2); % Mean MSE data for our Max Likelihood estimate

% 1(c)  Plotting MSE/ML with increasing measurements. As number of
% measurements increases, mean square error converges to 0.
figure;
plot(1:numMeas, BayesMSE,'b', 1:numMeas, ML_MSE,'r')
title('Bayes MMSE vs. ML: Convergence of MSE')
xlabel('Number of Measurements')
ylabel('Mean Square Error (MSE)')
ylim([0 10])
legend('Bayes MMSE', 'Max Likelihood')

% 1(d) - Plotting with theta pulled from different distributions
figure;
plot(1:numMeas, BayesMSE,'b', 1:numMeas, BayesMSE_uniform,'g', 1:numMeas,BayesMSE_exp, 'k')
title('Bayes MMSE: Correct vs Incorrect Prior')
xlabel('Number of Measurements')
ylabel('Mean Square Error (MSE)')
ylim([0 30])
legend('Gaussian (Correct)', 'Uniform (Incorrect)','Exponential (Incorrect)')

%% Scenario 2

% 2(a) - Mathematical model of BPSK signal w/ BPSK interference

% PDF of x, depending on point in time of signal, t, relative to start/end points of interference block, t1 t2
pdf = @(x,t,t1,t2,stdev) ... 
     (t<t1).*(normpdf(x,-1,stdev)/2 + normpdf(x,1,stdev)/2) + ... % BPSK signal only: 2 distributions, mean at +-1, variance is known
     (t>=t1).*(t<=t2).*(normpdf(x,-2,stdev)/4 + normpdf(x,0,stdev)/2+normpdf(x,2,stdev)/4) + ...  % This is when the interference is present. Mean can either be -1+-1, 1+1, or -1+1.
     (t>t2).*(normpdf(x,-1,stdev)/2 + normpdf(x,1,stdev)/2);  % Interference is not present. Signal only.

% Log likelihood fct = summation(ln(pdf))). ML picks values that maximize the log likelihood function.
lg = @(x,t,t1,t2,stdev) sum(log(pdf(x,t,t1,t2,stdev)));   


% 2(b) - Max Likelihood Estimator Implementation

% Set signal/noise parameters
t1 = 10;    % Interference is 1 contiguous block of noisefrom 
t2 = 80;    % Interference range: t1 to t2
nSymbols = 100;	% Number of symbols used in total
t = 1:nSymbols;   % Time

% Generate bits
signalBits = round(rand(1, nSymbols));     % Signal bits: ranges from symbols 1 to nSymbols
noiseBits = [zeros(1, t1) round(rand(1, t2 - t1)) zeros(1, nSymbols - t2)];  % Noise bits: ranges from t1 to t2, 0 elsewhere

% BPSK modulation: Bits must either be -1, 0, or 1
BPSKsignal = 2*(signalBits-0.5);    % BPSK modulate signal bits
BPSKnoise = 2*(noiseBits-0.5);      % BPSK modulate noise bits
BPSKnoise([1:t1*1, t2*1:end]) = 0;  % Revert bits outside interference range to 0
System = BPSKnoise + BPSKsignal;    % Add noise to signal

stDev = 0.01;   % Known standard deviation of noise (used to compute log likelihood)
x = System+stDev*randn(1,nSymbols);  % Generate one observation from signal+noise system (log likelihood of this will be computed)
curMaxLikelihood = -1e50;  % Very small starting value for likelihood, used for comparison.

for t1Sweep = 1:nSymbols
   for t2Sweep = t1Sweep:nSymbols
       loglikelihood = lg(x,t,t1Sweep,t2Sweep,stDev); % Compute log likelihood of x
       if loglikelihood > curMaxLikelihood     % Compare to current max likelihood
           curMaxLikelihood = loglikelihood;   % If larger, replace current max likelihood and update t1/t2 estimates
           t1Est = t1Sweep - 1;
           t2Est = t2Sweep + 1;
       end
   end
end

% 2(c) - Summary of Results

% Plot Signal with/without noise
figure;
plot(t, BPSKsignal);
legend('Pure BPSK Signal');
title('BPSK Modulated Signal');
xlabel('Symbols');
ylabel('Amplitude');
ylim([-2.5 2.5]);

figure;
plot(t, System);
legend('BPSK Signal + Interference');
title('BPSK Modulated Signal with Interference');
xlabel('Symbols');
ylabel('Amplitude');
ylim([-2.5 2.5]);

% Print ML estimates for parameters t1 and t2
t1Est 
t2Est

%% Scenario 3
