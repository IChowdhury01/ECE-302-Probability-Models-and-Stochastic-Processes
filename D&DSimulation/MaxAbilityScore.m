% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

clc, clear;
numIter = 100000;
AbilityScore = zeros(1,numIter);

for i = 1:numIter
    roll3d6 = randi(6,3,1); % Roll a 6 sided dice, 3 times. Record results in column vector
    AbilityScore(i) = sum(roll3d6);  % Sum up rolls to get 1 ability score
end

maxabilitycount = nnz(AbilityScore == 18); % Count number of times an ability score of 18 occurs
P_maxroll = maxabilitycount/numIter  