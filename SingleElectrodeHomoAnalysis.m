clearvars; clc;
load('InitialConditionsFull.mat')
load('solrbhomo.mat')
output = zeros(1,size(solrbhomo,1));

n =  0.5;
numactivated = zeros(1,size(solrbhomo.c,2));
solrbhomo.c = solrbhomo.c >= max(max(max(solrbhomo.c)))*n; % Must equal to at least half of all monte carle simulations 

for i = 1:size(solrbhomo.c,1) % For every step
    activated = squeeze(solrbhomo.c(i,:,:) > max(max(solrbhomo.c))*n);

end

output(i) = units(find(numactivated == 1,1));


