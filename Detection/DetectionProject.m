% Ivan Chowdhury, Josh Go
% ECE302 Probability Models & Stochastic Processes
% April 18, 2019
% Detection Project

clc;
clear all;
%% Question 2 - Introduction to pattern classification and machine learning
% Build a MAP classifier
% Divide data at random into 2 halves. One is for training, one is for
% testing.
% Training: run simulations to estimate mean/covariance
% Compute 4D pdf using mvnpdf. 
% Testing: Compute total Perror and confusion matrix


% Pre-Simulation
load('Iris.mat'); % Load Iris plant data set. 

nTrials = 100;   

% Preinitialize matrices
Perror = zeros(1,nTrials); % Probability of error per trial
ConfusionMatrix = zeros(3,3); % Confusion matrix

for i = 1:nTrials
% Split features into classes
    featuresC1 = features(1:50,:);
    featuresC2 = features(51:100,:);
    featuresC3 = features(101:150,:);
    
% Randomly Shuffle each class
    shuffleOrder = randperm(length(featuresC1));
    featuresC1 = featuresC1(shuffleOrder,:);
    featuresC2 = featuresC2(shuffleOrder,:);
    featuresC3 = featuresC3(shuffleOrder,:);

% Divide into training and testing data
    trainingFeaturesC1 = featuresC1(1:25,:);
    trainingFeaturesC2 = featuresC2(1:25,:);
    trainingFeaturesC3 = featuresC3(1:25,:);
    
    testFeatures = [featuresC1(26:50,:); featuresC2(26:50,:); featuresC3(26:50,:)];
    
% Training: estimate mean and covariance of training data
    mean1 = mean(trainingFeaturesC1);   % Mean vector of class 1 training features
    cov1 = cov(trainingFeaturesC1,1);   % Covariance matrix of class 1 training features
    
    mean2 = mean(trainingFeaturesC2);
    cov2 = cov(trainingFeaturesC2,1);
    
    mean3 = mean(trainingFeaturesC3);
    cov3 = cov(trainingFeaturesC3,1);
    
% Testing: Compute 4D PDF, Confusion matrix, and Probability of Error
    testCount = length(testFeatures);    % Number of test features used
    ErrorCount = 0;
    
    % Evaluate classifier on test features
    for j = 1:testCount
        curTest = testFeatures(j,:);  % Current test feature
        pdfC1 = mvnpdf(curTest,mean1,cov1); % Compute 4D pdf for each class of features
        pdfC2 = mvnpdf(curTest,mean2,cov2);
        pdfC3 = mvnpdf(curTest,mean3,cov3);
        
        % Compare PDFs. Choose the most likely classification for the test feature
        if (pdfC1/pdfC2 >= 1) && (pdfC1/pdfC3 >= 1) % Case 1: Classify test feature as class 1
            classifiedC1(j,:) = curTest;  
            
            % Verify whether classification was correct. Use result to compute confusion matrix
            if sum(ismember(features(1:50,:),curTest,'rows')) >= 1
                ConfusionMatrix(1,1) = ConfusionMatrix(1,1)+1; % Predicted Class 1 and its Class 1 (Detection)
            elseif sum(ismember(features(51:100,:),curTest,'rows')) >= 1
                ConfusionMatrix(2,1) = ConfusionMatrix(2,1)+1; % Predicted Class 1 but its Class 2 (Error)
                ErrorCount = ErrorCount + 1;
            else
                ConfusionMatrix(3,1) = ConfusionMatrix(3,1)+1; % Predicted Class 1 but its Class 3 (Error)
                ErrorCount = ErrorCount + 1;
            end
    
        elseif (pdfC2/pdfC1 >= 1) && (pdfC2/pdfC3 >= 1) % Case 2: Classify test as class 2
                classifiedC2(j,:) = curTest;
            if sum(ismember(features(1:50,:),curTest,'rows')) >= 1
                ConfusionMatrix(1,2) = ConfusionMatrix(1,2)+1; % Predicted Class 2 but its Class 1 (Error)
                ErrorCount = ErrorCount + 1;
            elseif sum(ismember(features(51:100,:),curTest,'rows')) >=1
                ConfusionMatrix(2,2) = ConfusionMatrix(2,2)+1; % Predicted Class 2 and its Class 2 (Detection)
            else
                ConfusionMatrix(3,2) = ConfusionMatrix(3,2)+1; % Predicted Class 2 but its Class 3 (Error)
                ErrorCount = ErrorCount + 1;
            end
        else    % Case 3: Classify as class 3.
            classifiedC3(j,:) = curTest;
             if sum(ismember(features(1:50,:),curTest,'rows')) >= 1
                ConfusionMatrix(1,3) = ConfusionMatrix(1,3)+1; % Predicted Class 3 but its Class 1 (Error)
                ErrorCount = ErrorCount + 1;
            elseif sum(ismember(features(51:100,:),curTest,'rows')) >= 1
                ConfusionMatrix(2,3) = ConfusionMatrix(2,3)+1; % Predicted Class 3 but its Class 2 (Error)
                ErrorCount = ErrorCount + 1;
            else
                ConfusionMatrix(3,3) = ConfusionMatrix(3,3)+1; % Predicted Class 3 and its Class 3 (Detection)
            end
        end
    end
    Perror(1,i) = ErrorCount/testCount;
end

% Final Results
ConfusionMatrix
PerrorTotal =  mean(Perror)