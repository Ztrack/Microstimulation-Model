clearvars; clc;

% THis code calculates the neuron rheobase, or the minimum amount of
% current required for a neuron to activate when taking into account only
% the center electrode.

% Features
calctype = 1; % If 1, this is the standard 4:1 ratio. if 2, this is a 5:X ratio with X being number of inhibitory
Motion_Axon = 0; % if set =1, enabled, connections between relevant motion neurons is established
Oscillatory_Behavior = 1; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location
axonswitch = 1; % if 0, then no axon is added. If 1, then axon is used in summation.

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
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
NumTrials = 100;
%inhibitoryfactor = [0.01 0.03 0.1 0.125 0.15 0.2]; % at rate = 40hz (for inhibitory), there is a X% inhibition factor active. This increases linearly with lambda.
inhibitoryfactor = [0.01];
NumInhibitoryMotif = 2; % Number of inhibitory neurons, if calctype = 2
neuron.lambdamod = lamdacombinedfun(neuron,zeros(size(neuron.type)),zeros(size(neuron.type)),1); % FR of neurons before any stimulation, includes synaptic connections

%% Loop Start
h = 50; % Steps
numrepeats = 3; % Number of overall repeats
ElectrodeDist = sqrt((sx/2-electrode.x).^2 + (sy/2-electrode.y).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode

for kk = 1:2
    lambdatype = kk;
    
    if lambdatype == 1
        %unitsmin = 1;
        unitsmax = 50;
    else
        %unitsmin = .0001;
        unitsmax = .001;
    end
    I0 = linspace(0,unitsmax,h);
    
    
    if calctype == 2
        load('InitialConditionsFullB.mat');
        NumNeuronsMotif = 5+NumInhibitoryMotif; % Reinitializes # Neurons per motif
        NumNeurons = NumMotifs*NumNeuronsMotif;
        a = zeros(1,NumInhibitoryMotif); b= ones(1,5-NumInhibitoryMotif); a = [a b]; a = [a 0 0 0 0 0];
        Extra_Inhib_Neurons = repmat(a,1,NumMotifs); Extra_Inhib_Neurons = find(Extra_Inhib_Neurons == 1); % Identifies extra neurons
        neuron.io.axon(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
        neuron.io.soma(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
        neuron.motif(Extra_Inhib_Neurons) = [];
        neuron.type(Extra_Inhib_Neurons) = [];
        neuron.oo.soma(Extra_Inhib_Neurons,:) = [];
        neuron.dirmult(Extra_Inhib_Neurons,:) = [];
        neuron.lambda = zeros(1,NumNeurons);
        neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
        neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
        neuron.inhibitory = find(neuron.type == 1); % New Inhibitory List
        neuron.excitatory = find(neuron.type == 2);
        neuron.neuron.lambdamod(neuron.type == 1) = 40;
        neuron.neuron.lambdamod(neuron.type == 2) = neuron.lambda(2) - neuron.lambda(1)*(neuron.lambda(1).*(inhibitoryfactor(1)/40));
        
    end
    
    Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
    
    parfor jj = 1:numrepeats
        
        Neuron_RB1 = NaN(1,NumNeurons);
        
        for ii = 1:length(I0)
            
            Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*neuron.dirmult(:,ElectrodeNo).*I0(ii) + neuron.io.axon(:,ElectrodeNo).*I0(ii).*axonswitch; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
            Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*I0(ii); % Summation of current directly from stimulus. AU irridance
            
            % Calculate neuron.lambda Change for poisson process
            [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,lambdatype);
            
            % Finding RB for each neuron
            
            for i = 1:NumNeurons
                if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
                    
                    Lambda_Hat_Spikes = Simple_PoissonGen(lambdahat(i), dt, NumTrials);
                    Y = prctile(Lambda_Hat_Spikes,bpct); % Calculates bottom xth percentile
                    if Y > neuron.lambdamod(i)+1
                        Neuron_RB1(i) = I0(ii);
                    end
                end
            end
        end
        Neuron_RB(jj,:) = Neuron_RB1;
    end
    
    %% Export
    if lambdatype == 1
        solrb.e1 = Neuron_RB;
        solrb.e.I0 = I0;
    elseif lambdatype == 2
        solrb.o1 = Neuron_RB;
        solrb.o.I0 = I0;
    end
    
end
if calctype == 1
    save('solrb1.mat','solrb','h','numrepeats','ElectrodeNo');
else
    save('solrb2.mat','solrb2','I0','h','numrepeats','ElectrodeNo');
end

