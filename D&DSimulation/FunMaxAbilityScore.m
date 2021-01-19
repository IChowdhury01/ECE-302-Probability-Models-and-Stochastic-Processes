% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

clc, clear;
numIter = 100000;
funAbilityScore = zeros(1, numIter);

for i = 1:numIter
    funRoll3d6 = sum(randi(6,3,3));
    funAbilityScore(i) = max(funRoll3d6);
end

maxAbilityCount = nnz(funAbilityScore == 18); % Count number of perfect ability scores in simulation.
P_funmaxroll = maxAbilityCount/numIter  