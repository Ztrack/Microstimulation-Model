clearvars; clc;
load('InitialConditionsFullB.mat')
% This initial condition has 5 Excitatory and 5 Inhibitory
IterationNumber = 0; % Number of inhibitory neurons

NumNeuronsMotif = 5+IterationNumber; % Reinitializes # Neurons per motif
NumNeurons = NumMotifs*NumNeuronsMotif;

a = zeros(1,IterationNumber); b= ones(1,5-IterationNumber); a = [a b]; a = [a 0 0 0 0 0];
Extra_Inhib_Neurons = repmat(a,1,NumMotifs); Extra_Inhib_Neurons = find(Extra_Inhib_Neurons == 1); % Identifies extra neurons
repmat(1:1:NumMotifs,NumNeuronsMotif,1); Neuron_Motif = Neuron_Motif(:)';
I0_Axon_Neurons(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
I0_Soma_Neurons(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
Neuron_Motif(Extra_Inhib_Neurons) = [];
Neuron_Type(Extra_Inhib_Neurons) = [];

%% Extracellular Current Matrix /Loop Start
h = 50; % Step size
I0 = h:h:20000;
numrepeats = 100;

Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes Lambda+1 spike
ElectrodeDist = sqrt((sx/2-ElectrodeX).^2 + (sy/2-ElectrodeY).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode

parfor jj = 1:numrepeats
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
for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if Neuron_Type(i) == 1
    Inhibitory_Effect(Neuron_Motif(i)) = Inhibitory_Effect(Neuron_Motif(i)) + Lambda_Hat(i).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
    Inhibitory_Effect_Lambda = FS_Lambda*Inhibitory_Factor*IterationNumber;
    end
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
NumTrials = 1000;
FS_Lambda_Spikes = mean(Simple_PoissonGen(Lambda(1), dt, NumTrials)); % Inhibitory lambda spikes
RS_Lambda_Spikes = mean(Simple_PoissonGen(Lambda(2), dt, NumTrials)); % Excitatory Lambda Spikes

RB_Exists1 = NaN(1,NumNeurons);
for i = 1:NumNeurons
    if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
        if Neuron_Type(i) == 1
            Lambda_Spikes = FS_Lambda_Spikes;
        else
            Lambda_Spikes = RS_Lambda_Spikes;
        end
        Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
        sort( Lambda_Hat_Spikes );
        Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
        if Y > mean(Lambda_Spikes+1)
            Neuron_RB1(i) = I0(ii);    
        end
    end
end
end
Neuron_RB(jj,:) = Neuron_RB1;
end

%% Export
if IterationNumber == 0
Neuron_RBB0 = Neuron_RB;
save('B0.mat','Neuron_RBB0');
elseif IterationNumber == 1
Neuron_RBB1 = Neuron_RB;
save('B1.mat','Neuron_RBB1');
elseif IterationNumber ==2
Neuron_RBB2 = Neuron_RB;
save('B2.mat','Neuron_RBB2');
elseif IterationNumber ==3
Neuron_RBB3 = Neuron_RB;
save('B3.mat','Neuron_RBB3');
elseif IterationNumber ==4
Neuron_RBB4 = Neuron_RB;
save('B4.mat','Neuron_RBB4');
elseif IterationNumber ==5
Neuron_RBB5 = Neuron_RB;
save('B5.mat','Neuron_RBB5');
end

%% Functions

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random(:,:) < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end