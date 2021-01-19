% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

clc, clear;
numIter = 100000;
numTrolls = 6;

trollPack = zeros(1,numTrolls);
remainingHP = zeros(numIter,numTrolls);
numAlive = zeros(numIter,1);
survivorHP = zeros(1,numIter);

for i = 1:numIter
    for j = 1:6
        trollPack(j) = randi(4,1,1);    % Randomly generate troll pack
    end
    
    FIREBALL = sum(randi(2,2,1)); % Randomly generate fireball
    
    remainingHP(i,:) = trollPack - FIREBALL; % Find remaining troll HP
    numAlive(i) = nnz(remainingHP(i,:) >= 1); % Counts how many trolls from each pack is still alive
    if (numAlive(i) == 1)
        survivorIndex = find(remainingHP(i,:) >= 1); % Find HP of survivor troll
        survivorHP(i) = remainingHP(i,survivorIndex); % skipped indices in matrix will have zero values, these will be ignored
    end
end

expectedSurvivorHP = mean(nonzeros(survivorHP)) 