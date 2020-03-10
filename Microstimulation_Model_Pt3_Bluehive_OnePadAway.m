clearvars; clc;
load('InitialConditionsFull.mat')



%% Extracellular Current Matrix /Loop Start

Neuron_Onepadaway = find(Pad_Neuron == 3 | Pad_Neuron == 7 | Pad_Neuron == 9 | Pad_Neuron == 13);

%h = 50; % Step size
I0 = 880; % CH for pad D3
numrepeats = 5000;

Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes Lambda+1 spike
ElectrodeDist = sqrt((sx/2-ElectrodeX).^2 + (sy/2-ElectrodeY).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode

for jj = 1:numrepeats
Neuron_RB1 = NaN(1,NumNeurons);
Ie_Axon_Neurons = zeros(1,length(ElectrodeNo));
Ie_Soma_Neurons = zeros(1,length(ElectrodeNo));

for ii = 1:length(I0)

Ie_Axon_Neurons = sum(I0_Axon_Neurons(:,ElectrodeNo).*I0(ii),2);
Ie_Soma_Neurons = sum(I0_Soma_Neurons(:,ElectrodeNo).*I0(ii),2);
Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

%% Multiplier Matrices
dt = 1/1000; % step/bin time length 1ms
t_total= 0:dt:T-dt ;
NumSteps = 1/dt; % Number of steps/bins in time t_total
FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection
Lambda = zeros(NumNeurons,1);

for i = 1:NumNeurons
    if Neuron_Type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
        Lambda(i) = FS_Lambda; % Lambda regular, no microstimulation
    end
end

Inhibitory_Effect = zeros(NumMotifs,1);
Inhibitory_Effect_Lambda = zeros(NumMotifs,1);
for i = 1:NumMotifs
    Inhibitory_Effect(i) = Lambda_Hat(Neuron_Inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
    Inhibitory_Effect_Lambda = FS_Lambda*Inhibitory_Factor;
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if Neuron_Type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(Neuron_Motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
        Lambda(i) = RS_Lambda - Inhibitory_Effect_Lambda; % Lambda regular, no microstimulation
    else
        % Inhib Neurons do not inhib themselves
    end
end

%% Finding RB for each neuron
Lambda_Spikes = [0 0];
Lambda_Spikes(1) = FS_Lambda; % Inhibitory lambda spikesLambda_Spikes
Lambda_Spikes(2) = RS_Lambda; % Excitatory Lambda Spikes
NumTrials = 1000;

for i = 1:length(Neuron_Onepadaway)
    if isnan(Neuron_RB1(Neuron_Onepadaway(i))) % If RB does not exist, continue, otherwise this neuron is skipped
%         if sum(Neuron_Inhibitory_Oscillatory == i) > 0 % If the neuron has oscillatory behavior then use:
%             Lambda_Hat_Spikes = Oscillatory_PoissonGen(Lambda_Hat(i), dt, NumTrials);
%         else % Otherwise use the simple function:
            Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(Neuron_Onepadaway(i)), dt, NumTrials);
%         end
        
        sort( Lambda_Hat_Spikes );
        Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
        if Y > mean(Lambda_Spikes(Neuron_Type(Neuron_Onepadaway(i)))+1)
            Neuron_RB1(Neuron_Onepadaway(i)) = 1;    
        end
    end
end
end
Neuron_RB(jj,:) = Neuron_RB1;
end

%% Export

Neuron_RB_Onepadaway = Neuron_RB;
save('Onepadaway.mat','Neuron_RB_Onepadaway');

%% Functions

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random(:,:) < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end

function Trial_Spikes = Oscillatory_PoissonGen(freqOscillation, dt, NumTrials)

NumSteps = 1/dt;
Trial_Spikes = nan(NumTrials,1);
t = 1:1:1000;

for n = 1 : NumTrials
    
    signal = 1 + (0-1).*rand(NumSteps,1);
    
    curFreq = randi([1, length(freqOscillation)],1,1);
    oscillationModulator = sin(2 * pi * freqOscillation(curFreq) * t/1000);
    oscillationModulator(oscillationModulator < 0) = 0;
    
    %-------
    %adding a bit of normal dist noise to the underlying oscillation function
    x = normrnd(0,5,[length(oscillationModulator), 1]);
    x = x - min(x);
    x = (x)./max(x(:));
    oscillationModulator = oscillationModulator .* x';
    oscillationModulator = oscillationModulator ./ max(oscillationModulator(:));
    modSignal = oscillationModulator';
    
    spikes = (signal < (modSignal/4.5));
    Trial_Spikes(n) = sum(spikes);
end
end