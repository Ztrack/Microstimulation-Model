clearvars; clc;

% THis code calculates the neuron rheobase, or the minimum amount of
% current required for a neuron to activate when taking into account only
% the center electrode.

% Features
calctype = 1;
NumInhibitoryMotif = 2; % Number of inhibitory neurons, if calctype = 2
Oscillatory_Behavior = 0; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location
axonswitch = 1; % if 0, then no axon is added. If 1, then axon is used in summation.

% Apply Features
load('InitialConditionsFull.mat') % 4:1 Excitatory to inhibitory ratio
if Directional_Current_Modifier == 0 % Directionality component to current summation (Based on axon hilic)
    neuron.dirmult(:,:) = 1;
end

% Parameters
h = 50; % Steps
numrepeats = 3; % Number of overall repeats
NumTrials = 20;
bpct = 05; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
neuron.lambda = zeros(params.numneurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% Calculations
%neuron.lambdamod = lamdacombinedfun(neuron,params,zeros(size(neuron.type)),zeros(size(neuron.type)),1); % FR of neurons before any stimulation, includes synaptic connections
parfor i = 1:NumTrials
    lambdamod(i,:) = IntegrateAndFire(neuron,params,zeros(params.numneurons,1)); % FR of neurons before any stimulation, includes synaptic connections
end
neuron.lambdamod = mean(lambdamod,1);

%% Loop Start
solrb2.e1 = NaN(length(electrode.x),numrepeats,params.numneurons); % Current RB init
solrb2.o1 = NaN(length(electrode.x),numrepeats,params.numneurons); % Opto RB init

for kk = 1:length(electrode.x)
    for jj = 1:numrepeats
        
        for lambdatype = 1:2
            
            Neuron_RB = NaN(1,params.numneurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
            
            if lambdatype == 1
                unitsmin = 0;
                unitsmax = 1200;
            else
                unitsmin = 0;
                unitsmax = 0.25;
            end
            I0 = linspace(unitsmin,unitsmax,h);
            
            for ii = 1:length(I0)
                
                Ie_Neurons = neuron.io.soma(:,kk).*neuron.dirmult(:,kk).*I0(ii) + neuron.io.axon(:,kk).*I0(ii).*axonswitch; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
                Il_Neurons = neuron.oo.soma(:,kk).*I0(ii); % Summation of current directly from stimulus. AU irridance
                
                % Calculate neuron.lambda Change
                [lambdahat] = lamdacombinedfun(neuron,params,Ie_Neurons,Il_Neurons,lambdatype); % Change in firing rate (Hz) due to electrical stimulation or optogenetics
                
                % Run spike generator
                out_FR = zeros(NumTrials,params.numneurons);
                parfor i = 1:NumTrials
                    out_FR(i,:) = IntegrateAndFire(neuron,params,lambdahat);
                end
                
                % Determine if neuron spiked 1 above baseline
                for i = 1:params.numneurons
                    if isnan(Neuron_RB(i)) % Continue only if nan, saves calculation time
                        Y = prctile(out_FR(:,i),bpct); % Calculates bottom xth percentile
                        if Y > neuron.lambdamod(i)+1
                            Neuron_RB(i) = I0(ii); % Update Neuron Rhoebase
                        end
                    end
                end
                if sum(isnan(Neuron_RB)) <= params.numneurons*.20
                    break
                end
            end
            
            % Save Data
            if lambdatype == 1
                solrb2.e1(kk,jj,:) = Neuron_RB;
                solrb2.e.I0 = I0;
                
            elseif lambdatype == 2
                solrb2.o1(kk,jj,:) = Neuron_RB;
                solrb2.o.I0 = I0;
            end
        end
        
    end
end

%% Export
save('solrb2.mat','solrb2');
