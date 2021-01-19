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

for i = 1:numIter
    for j = 1:6
        trollPack(j) = randi(4,1,1);    % Randomly generate troll pack
    end
    
    FIREBALL = sum(randi(2,2,1)); % Randomly generate fireball

    remainingHP(i,:) = trollPack - FIREBALL; % Find remaining troll HP
    numAlive(i,1) = nnz(remainingHP(i,:) > 0); % Counts how many trolls from each pack is still alive
end

numDeadPacks = nnz(numAlive == 0); % Count how many troll packs were fully killed
P_slayall = numDeadPacks/numIter % P(all trolls slayed)