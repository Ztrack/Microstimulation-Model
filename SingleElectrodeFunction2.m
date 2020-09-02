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
simulation = 4; % Number of overall repeats for monte carle simulation
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

%% Loop Start
h = 50; % number of steps for current/LI

% if lambdatype == 1
%     unitsmax = 30000; % Point at which 100% of neurons are activated
%     units50 = 10000; % Point at which 50% of neurons are activated
% else
%     unitsmax = 100000;
%     units50 = 30000; 
% end
% units = linspace(0,units50,h*.8);  % Current OR liminous intensity Steps
% units = [units linspace(units50+unitsmax*.2*h,unitsmax,h*.2)];

if lambdatype == 1
    unitsmin = 1;
    unitsmax = 50; 
else
    unitsmin = .0001;
    unitsmax = .05; 
end
units = linspace(0,unitsmax,h); 

output = zeros(length(electrode.x),h,NumNeurons); % Stores poisson spike rate for every output, logical matrix

parfor j = 1:length(electrode.x) % Iterate every electrode
    output1 = zeros(h,NumNeurons);
    ElectrodeNo = j;
    
    for i = 1:length(units)
        
        Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*neuron.dirmult(:,ElectrodeNo).*units(i) + neuron.io.axon(:,ElectrodeNo).*units(i); % Summation of current directly from stimulus + backpropogated up by axons. AU Current
        Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*units(i); % Summation of current directly from stimulus. AU luminous intensity
        
        % Calculate neuron.lambda Change for poisson process
        [lambdahat,lambdamod] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,lambdatype);
        BestNeurons = lambdahat >= median(lambdahat); % Calculate only the easiest to activate 50%
        
        % Calculate Poisson spike rates for each neuron
        for ii = 1:NumNeurons
            if BestNeurons(ii) == 1
                y = Simple_PoissonGen2(lambdahat(ii), dt, NumTrials,simulation,bpct) > lambdamod(ii)+1; % Calculate Lambda
                output1(i,ii) = sum(y);
                %outputactivated(ii) = sum(y) >= simulation.*(.5); % Useful for debugging
            end
        end
        
    end
    
    output(j,:,:) = output1; % Output result for this parfor loop
    disp(['Electrode ' num2str(j) ' done']);
end

if lambdatype == 1
    save('SEoutputc.mat','output','units','-v7.3'); % Single electrode results output
elseif lambdatype == 2
    save('SEoutputo.mat','output','units','-v7.3'); % Single electrode results output
end
