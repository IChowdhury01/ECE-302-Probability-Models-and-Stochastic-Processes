% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

clc, clear;
numIter = 100000;

trollHP = zeros(1,numIter);
FIREBALL = zeros(1,numIter);

for i = 1:numIter
    trollHP(i) = randi(4,1,1); % Roll a 4 sided dice, 1 time
    FIREBALL(i) = sum(randi(2,2,1)); % Roll a 2-sided dice, twice, and sum results.
end

% Averages
avgHP = mean(trollHP)
avgFIREBALLdmg = mean(FIREBALL)

% Calculate PMF
P_trollHP1 = nnz(trollHP == 1)/numIter % P(Troll HP = 1)
P_trollHP2 = nnz(trollHP == 2)/numIter % P(Troll HP = 2)
P_trollHP3 = nnz(trollHP == 3)/numIter % P(Troll HP = 3)
P_trollHP4 = nnz(trollHP == 4)/numIter % P(Troll HP = 4)

P_fireball2 = nnz(FIREBALL == 2)/numIter % P(FIREBALL DMG = 2)
P_fireball3 = nnz(FIREBALL == 3)/numIter % P(FIREBALL DMG = 3)
P_fireball4 = nnz(FIREBALL == 4)/numIter % P(FIREBALL DMG = 4)

HP = [1 2 3 4];
FireBallDMG = [2 3 4];
P_TrollHP = [P_trollHP1 P_trollHP2 P_trollHP3 P_trollHP4];
P_fireball = [P_fireball2 P_fireball3 P_fireball4];

figure;
bar(HP,P_TrollHP);
title("Troll HP PMF");
xlabel("Troll HP");
ylabel("Probability");

figure;
bar(FireBallDMG,P_fireball);
title("Fireball DMG PMF");
xlabel("Fireball DMG");
ylabel("Probability");
