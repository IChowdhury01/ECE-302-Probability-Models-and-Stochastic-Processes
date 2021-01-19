% Ivan Chowdhury
% Stochastics
% Spring 2019
% Dungeons and Dragons Simulations

clc, clear;
numIter = 100000;
totalDMG = zeros(1, numIter);

for i = 1:numIter
    attackRoll = randi(20);
    if (attackRoll >= 11)
       SwordDMG = sum(randi(6,1,2));
       totalDMG(i) = SwordDMG;
       
       attackRoll2 = randi(20);
       if (attackRoll2 >= 11)
          HammerDMG = randi(4);
          totalDMG(i) = SwordDMG + HammerDMG;
       end
    end
end

expectedDMG = mean(totalDMG)