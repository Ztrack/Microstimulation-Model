clearvars; clc;
load('InitialConditionsFull.mat')
load('solrbhomo.mat')
output = zeros(1,size(solrbhomo,1));

n =  0.5;
% solrbhomo.c = Current/opto steps x number of neurons x time period
solrbhomo.c = solrbhomo.c >= max(max(max(solrbhomo.c)))*n; % Must equal to at least half of all monte carle simulations 
Neuron_RB = nan(NumNeurons,10); % 10 Rheobase values for every neuron, representing 10 sections of time

for ii = 1:size(solrbhomo.c,2) % For every neuron
    for iii = 1:size(solrbhomo.c,3) % For every time period
        
        RB = units(find(solrbhomo.c(:,ii,iii),1)); % Find the first instance that this neuron activated at this time period
        
        if RB > 0 % If there is a rhoebase for this combination, record it
            Neuron_RB = RB; % 
        else
            Neuron_RB = NaN; % Otherwise, there is no solution
        end
    end
end

% Now that we know the Rhoebase for every neuron & time period, we create a
% vector y which stores the percent of all neurons activated across units
% vector.

for i = 1:length(solrbhomo.units_c)
    stepsol.current.all(:,i) = sum(Neuron_RB                          <=solrbhomo.units_c(i),2)*100/NumNeurons;
    stepsol.current.excitatory(:,i) = sum(Neuron_RB(neuron.excitatory)     <=solrbhomo.units_c(i),2)*100/length(neuron.excitatory);
    stepsol.current.inhibitory(:,i) = sum(Neuron_RB(neuron.inhibitory)     <=solrbhomo.units_c(i),2)*100/length(neuron.inhibitory);
    stepsol.current.motion(:,i) = sum(Neuron_RB(neuron.motion.number)    <=solrbhomo.units_c(i),2)*100/length(neuron.motion.number);
    stepsol.current.nonmotion(:,i) = sum(Neuron_RB(neuron.nonmotion.number) <=solrbhomo.units_c(i),2)*100/length(neuron.nonmotion.number);
end

