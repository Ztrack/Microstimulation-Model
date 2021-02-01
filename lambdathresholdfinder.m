function [lambdathreshold] = lambdathresholdfinder(neuron,CI)

% This function finds the threshold value at which a neuron will fire once
% over the baseline given a confidence interval CI. This is done to save
% computational time during the optimization function. 
% The function tests all values from lambda1 to lambda2 using the poisson
% spike train generator using the res resolution parameter. Typical value
% for res is rounded to the nearest 0.01. 

NumTrials = 20;
dt = 1/1000;
lambda = round(min(neuron.lambdamod),1):0.1:round(max(neuron.lambdamod*2),1);
Y = nan(size(lambda));
lambdathreshold = nan(length(neuron.lambdamod),1);

parfor i = 1:length(lambda)
    
    Lambda_Hat_Spikes = Simple_PoissonGen(lambda(i), dt, NumTrials); % Generates NumTrials number of poisson spike trains
    Y(i) = prctile(Lambda_Hat_Spikes,100-CI); % Calculates bottom xth percentile & stores output

end

% Now that the lookup table is complete, we connect the corresponding
% threshold for every neuron.

for i = 1:length(neuron.lambdamod)
    
    lambdathreshold(i) = lambda(find(Y > neuron.lambdamod(i)+1,1)); % Find the lambda value that results in a Y greater than lambda+1
    
end

end