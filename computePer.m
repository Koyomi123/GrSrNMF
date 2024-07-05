function [CA, CR, CP, CF, CNMI] = computePer(Z,GT_Cube)
T=length(GT_Cube);
for t=1:T
    gt{1}=GT_Cube{t};
    [CA(t), CR(t), CP(t), CF(t), CNMI(t)] = computePerformance(Z(:,t),gt);
end