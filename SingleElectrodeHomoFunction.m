clearvars; clc;

% THis code calculates the neuron rheobase, or the minimum amount of
% current required for a neuron to activate when taking into account only
% the center electrode.

% Features
calctype = 1; % If 1, this is the standard 4:1 ratio. if 2, this is a 5:X ratio with X being number of inhibitory
Motion_Axon = 0; % if set =1, enabled, connections between relevant motion neurons is established
Oscillatory_Behavior = 1; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location
lambdatype = 1; % What type of calcultion stimulus is presented. 1= current, 2 = opsin

% Apply Features
if calctype == 1 %load initial condition
    load('InitialConditionsFull.mat') % 4:1 Excitatory to inhibitory ratio
else
    load('InitialConditionsFullB.mat') % Multiple Inhibitory neurons present
end
if Directional_Current_Modifier == 0 % Directionality component to current summation (Based on axon hilic)
    neuron.dirmult(:,:) = 1;
end
if Motion_Axon == 0
    I0_Motion_Neurons = zeros(NumNeurons,length(electrode.x));
end

% Parameters
bpct = 05; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
NumTrials = 50; % Number of trials per poisson process
simulation = 20; % Number of overall repeats for monte carle simulation
inhibitoryfactor = 0.01; % at rate = 40hz (for inhibitory), there is a X% inhibition factor active. This increases linearly with lambda.
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
[~,lambdamod] = lamdacombinedfun(neuron,zeros(1,1000),zeros(1,1000),1);
%% Lambda Homo Calc
% Calculating the number of spikes in each region of time, to later analyze
% if we reached +1 spikes within this region.


neuron.adapt.lambdahomo = zeros(3,NumNeurons,10);
for i = 1:3
    for ii = 1:NumNeurons

        neuron.adapt.lambdahomo(i,ii,:) = mean(nonhomoPoissonGen(lambdamod(ii), dt, NumTrials,simulation,bpct,neuron.adapt.ratefunction(i,:)),1); % Sum of these numbers = input lambda, or close to

    end
end

%% Loop Start
h = 200; % number of steps for current/LI

% if lambdatype == 1
%     unitsmax = 30000; % Point at which 100% of neurons are activated
%     units50 = 10000; % Point at which 50% of neurons are activated
% else
%     unitsmax = 100000;
%     units50 = 30000;
% end
% units = linspace(0,units50,h*.8);  % Current OR liminous intensity Steps
% units = [units linspace(units50+unitsmax*.2*h,unitsmax,h*.2)];


ElectrodeNo = 45; % Number of the electrode, 45 = center
numrepeats = 3;

for j = 1:2
    
    lambdatype = j;
    if lambdatype == 1
        unitsmin = 1;
        unitsmax = 50;
    else
        unitsmin = .0001;
        unitsmax = .04;
    end
    units = linspace(0,unitsmax,h);
    output = zeros(h,NumNeurons,numrepeats,10); % Stores poisson spike rate for every output, logical matrix
    
    parfor i = 1:length(units)
        
        output1 = zeros(NumNeurons,numrepeats,10);
        Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*neuron.dirmult(:,ElectrodeNo).*units(i) + neuron.io.axon(:,ElectrodeNo).*units(i); % Summation of current directly from stimulus + backpropogated up by axons. AU Current
        Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*units(i); % Summation of current directly from stimulus. AU luminous intensity
        
        % Calculate neuron.lambda Change for poisson process
        [lambdahat,lambdamod] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,lambdatype);
        
        % Calculate Poisson spike rates for each neuron
        for ii = 1:NumNeurons
            for iii = 1:numrepeats % Repetition
                x = neuron.adapt.lambdahomo(neuron.adapt.type(ii),ii,:)+1;
                x = reshape(x,1,10);
                y = nonhomoPoissonGen(lambdahat(ii), dt, NumTrials,simulation,bpct,neuron.adapt.ratefunction(neuron.adapt.type(ii),:)); % Calculate Lambda
                output1(ii,iii,:) = sum(y > x);
                %outputactivated(ii) = sum(y) >= simulation.*(.5); % Useful for debugging
            end
        end
        
        output(i,:,:,:) = output1;
        
    end
    
    if lambdatype == 1
        solrbhomo.c = output;
        solrbhomo.units_c = units;
    elseif lambdatype == 2
        solrbhomo.o = output;
        solrbhomo.units_o = units;
    end
    
end
save('solrbhomo.mat','solrbhomo','-v7.3'); % Single electrode results output
