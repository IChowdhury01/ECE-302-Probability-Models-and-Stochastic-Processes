% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

% Very high number of iterations needed to generate a perfect char

clc, clear;
numIter = 100000;
charStats = zeros(numIter,6);
FontaineCtr = 0;

for i = 1:numIter
    for j = 1:6
        funRoll3d6 = sum(randi(6,3));
        charStats(i,j) = max(funRoll3d6); % Assign 6 ability scores per character (row) 
    end
    
    numStat9(i) = nnz(charStats(i,:)==18); % Counts how many of each character's stats is equal to 18
    if numStat9(i) == 6
        FontaineCtr = FontaineCtr + 1;
    end
end

P_perfectChar = FontaineCtr/numIter  
