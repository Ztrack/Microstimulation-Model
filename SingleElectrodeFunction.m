clearvars; clc;

% THis code calculates the neuron rheobase, or the minimum amount of
% current required for a neuron to activate when taking into account only
% the center electrode.

% Features
axonswitch = 1; % if 0, then no axon is added. If 1, then axon is used in summation.
neuronweights = 1; % If neuron weights apply, leave as 1.
inhibitoryweights = 1; % If neuron inhibitory weights apply, leave as 1.

% Apply Features
load('InitialConditionsFull.mat')
if neuronweights == 0
    neuron.weight.matrix = zeros(params.numneurons,params.numneurons);
end
if inhibitoryweights == 0
    neuron.weight.matrix(neuron.weight.matrix==neuron.weight.WEI) = 0; % Inhibitory onto excitatory connections will = 0
end

% Parameters
h = 48; % Steps
numrepeats = 1; % Number of overall repeats
NumTrials = 50;
bpct = 05; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
neuron.lambda = zeros(params.numneurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% Calculations

parfor i = 1:NumTrials
    lambdamod(i,:) = IntegrateAndFire(neuron,params,zeros(params.numneurons,1)); % FR of neurons before any stimulation, includes synaptic connections
end
neuron.lambdamod = mean(lambdamod,1);
%ElectrodeDist = sqrt((params.sx/2-electrode.x).^2 + (params.sy/2-electrode.y).^2);
%ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode
ElectrodeNo = 45;

%% Loop Start
solrb.e1 = NaN(numrepeats,params.numneurons); % Current RB init
solrb.o1 = NaN(numrepeats,params.numneurons); % Opto RB init

for jj = 1:numrepeats
    
    for lambdatype = 1:2
        Neuron_RB = NaN(1,params.numneurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
        
        if lambdatype == 1
            unitsmin = 0;
            unitsmax = 20;
        else
            unitsmin = 0;
            unitsmax = 0.015;
        end
        I0 = linspace(unitsmin,unitsmax,h);
        
        for ii = 1:length(I0)
            
            Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*I0(ii) + neuron.io.axon(:,ElectrodeNo).*I0(ii).*axonswitch; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
            Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*I0(ii); % Summation of current directly from stimulus. AU irridance
            
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
        end
        
        % Save Data
        if lambdatype == 1
            solrb.e1(jj,:) = Neuron_RB;
            solrb.e.I0 = I0;
            
        elseif lambdatype == 2
            solrb.o1(jj,:) = Neuron_RB;
            solrb.o.I0 = I0;
        end
    end
    
end


%% Export
save('solrb1.mat','solrb','numrepeats','ElectrodeNo');
