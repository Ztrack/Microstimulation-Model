clearvars; clc;

% Features
calctype = 1;
Motion_Axon = 0; % if set =1, enabled, connections between relevant motion neurons is established
Oscillatory_Behavior = 1; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location

% Apply Features
if calctype == 1 %load initial condition
    load('InitialConditionsFull.mat')
else
    load('InitialConditionsFullB.mat')
end
if Directional_Current_Modifier == 0
    neuron.dirmult(:,:) = 1;
end
if Motion_Axon == 0
    I0_Motion_Neurons = zeros(NumNeurons,length(electrode.x));
end

% Parameters
Axonal_Mult = .01;
bpct = 50; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
%neuron.inhibitoryfactor = [0.01 0.03 0.1 0.125 0.15 0.2];
neuron.inhibitoryfactor = [0.01];
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
NumTrials = 100;

%% Loop Start
h = 500; % Step size
I0 = 0:h:100000;  % Current Steps
numrepeats = 12; % Number of overall repeats
ElectrodeDist = sqrt((sx/2-electrode.x).^2 + (sy/2-electrode.y).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode


for kk = 1:length(neuron.inhibitoryfactor)
    
    Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
    
    parfor jj = 1:numrepeats
        
        Neuron_RB1 = NaN(1,NumNeurons);
        
        for ii = 1:length(I0)
            Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
            Ie_Soma_Neurons = zeros(NumNeurons,1);
            
            for i = 1:length(ElectrodeNo) % Summation of current for every neuron component by every electrode & its corresponding current
                Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*I0(ii);
                Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*I0(ii).*neuron.dirmult(:,ElectrodeNo);
            end
            
            Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons
            
            % Calculate neuron.lambda Change for poisson process
            
            lambdahat = lamdafun(neuron,Ie_Soma_Axon_Neurons,kk);
            
            % Finding RB for each neuron
            
            for i = 1:NumNeurons
                if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
                    if neuron.oscillatorytype(i) == 1 & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
                        Lambda_Hat_Spikes = Oscillatory_PoissonGen(lambdahat(i), dt, NumTrials);
                    else % Otherwise use the simple function:
                        Lambda_Hat_Spikes = Simple_PoissonGen(lambdahat(i), dt, NumTrials);
                    end
                    
                    Y = prctile(Lambda_Hat_Spikes,bpct); % Calculates bottom xth percentile
                    if Y > mean(neuron.lambda(i)+1)
                        Neuron_RB1(i) = I0(ii);
                    end
                end
            end
        end
        Neuron_RB(jj,:) = Neuron_RB1;
    end
    
    %% Export
    if kk == 1
        solrb.a1 = Neuron_RB;
    elseif kk ==2
        solrb.a2 = Neuron_RB;
    elseif kk ==3
        solrb.a3 = Neuron_RB;
    elseif kk ==4
        solrb.a4 = Neuron_RB;
    elseif kk ==5
        solrb.a5 = Neuron_RB;
    elseif kk ==6
        solrb.a5 = Neuron_RB;
    end
    
end

if calctype == 1
    save('solrb1.mat','solrb','I0','h','numrepeats','ElectrodeNo');
else
    save('solrb2.mat','solrb','I0','h','numrepeats','ElectrodeNo');
end

