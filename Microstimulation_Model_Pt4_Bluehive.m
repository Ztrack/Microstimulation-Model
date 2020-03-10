clearvars; clc;
load('InitialConditionsFull.mat')

%%
% For this first part, we only care about the motion neurons we are trying
% to stimulate. We want to find their RB values for use in the next part.

I0_Total = I0_Soma_Neurons + I0_Axon_Neurons;
a = I0_Total(:);
a = sort(a,'descend');
for i = 1:45
    b(i) = find(I0_Total == a(i));
end
[Neuron_Target(1,:),Neuron_Target(2,:)] = ind2sub(size(I0_Soma_Neurons), b);

%% Finding the best electrode-neuron pairs
h = 1; % Step size
I0 = 1:h:10;
numrepeats = 2;
Inhibitory_Factor = 0.0; % Control
ElectrodeNo = 1:100; % Activate all electrodes, we want to activate every motion neuron
Neuron_Motion_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes Lambda+1 spike

for jj = 1:numrepeats
Neuron_RB1 = NaN(1,NumNeurons);
Ie_Axon_Neurons = zeros(1,length(ElectrodeNo));
Ie_Soma_Neurons = zeros(1,length(ElectrodeNo));

for ii = 1:length(I0)

Ie_Axon_Neurons = sum(I0_Axon_Neurons(:,ElectrodeNo).*I0(ii),2);
Ie_Soma_Neurons = sum(I0_Soma_Neurons(:,ElectrodeNo).*I0(ii),2);

Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

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

Lambda_Spikes = [0 0];
Lambda_Spikes(1) = FS_Lambda; % Inhibitory lambda spikesLambda_Spikes
Lambda_Spikes(2) = RS_Lambda; % Excitatory Lambda Spikes
NumTrials = 1000;

for i = 1:length(Neuron_Motion)
    if isnan(Neuron_RB1(Neuron_Motion(1,i))) % If RB does not exist, continue, otherwise this neuron is skipped
%         if sum(Neuron_Inhibitory_Oscillatory == i) > 0 % If the neuron has oscillatory behavior then use:
%             Lambda_Hat_Spikes = Oscillatory_PoissonGen(Lambda_Hat(i), dt, NumTrials);
%         else % Otherwise use the simple function:
            Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(Neuron_Motion(1,i)), dt, NumTrials);
%         end
        
        sort( Lambda_Hat_Spikes );
        Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
        if Y > mean(Lambda_Spikes(Neuron_Type(Neuron_Motion(1,end)))+1)
            Neuron_RB1(Neuron_Motion(1,i)) = I0(ii);    
        end
    end
end
end
Neuron_Motion_RB(jj,:) = Neuron_RB1;
end

NeuronStimDist = NaN(NumNeurons,length(ElectrodeX));
Neuron_Electrode = zeros(NumNeurons,3);
for i = 1:length(ElectrodeX)
    for j = 1:NumNeurons
        NeuronStimDist(j,i) = sqrt((NeuronX(j)-ElectrodeX(i)).^2 + (NeuronY(j) - ElectrodeY(i)).^2);
    end
end
for i = 1:NumNeurons
    a = min(NeuronStimDist(i,:));
    Neuron_Electrode(i,1) = i; % Neuron Number
    Neuron_Electrode(i,2) = find(NeuronStimDist(i,:) == a,1); % Closest electrode to this neuron. In rare case, 1 or 2 neurons are equidistant to 2 electrodes
    Neuron_Electrode(i,3) = a; % Stores length to electrode
    Neuron_Electrode(i,4) = nanmean(Neuron_Motion_RB(:,i)); % Neuron RB Value
end
Neuron_Motion_Electrode = Neuron_Electrode(Neuron_Motion(1,:),:);
a = []; k = 1;
for i = 1:NumNeurons % NaN Filter - remove any neurons not activated
    if isnan(Neuron_Electrode(i,4))
        %
    else
        a(k,:) = Neuron_Electrode(i,:);
        k = k+1;
    end
end
        
Neuron_Electrode_Targeted = sortrows(a,4,'ascend');
Neuron_Electrode_Targeted = Neuron_Electrode_Targeted(1:45,:); % We target the first 45 motion tuned neurons
Stim_Electrode = [];
Stim_Electrode(:,1) = Neuron_Electrode_Targeted(:,2); % Electrode # , Current intensity
Stim_Electrode(:,2) = Neuron_Electrode_Targeted(:,4).*2; % Current intensity * 2 = Chronaxi
for i = 1:length(Stim_Electrode)
    if sum(Stim_Electrode(i,1) == Stim_Electrode(:,1)) > 1 
    Stim_Electrode(i,2) = max(Stim_Electrode(Stim_Electrode(i,1) == Stim_Electrode(:,1),2));
    end
end
Stim_Electrode = unique(sort(Stim_Electrode,1), 'rows'); % Removes repeated electrode numbers

%% Determine which other neurons were activated

numrepeats = 1;
Neuron_Activated = zeros(1,NumNeurons); % 0 = inactivated, 1 = activated
Ie_Axon_Neurons = zeros(NumNeurons,1);
Ie_Soma_Neurons = zeros(NumNeurons,1);


for i = 1:length(Stim_Electrode) % Summation of current for every neuron component by every electrode & its corresponding current
Ie_Axon_Neurons = Ie_Axon_Neurons + I0_Axon_Neurons(:,Stim_Electrode(i,1)).*Stim_Electrode(i,2);
Ie_Soma_Neurons = Ie_Soma_Neurons + I0_Soma_Neurons(:,Stim_Electrode(i,1)).*Stim_Electrode(i,2);
end

Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

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

Lambda_Spikes = [0 0];
Lambda_Spikes(1) = FS_Lambda; % Inhibitory lambda spikesLambda_Spikes
Lambda_Spikes(2) = RS_Lambda; % Excitatory Lambda Spikes
NumTrials = 5000;


%Neuron_Activated1 = zeros(1,NumNeurons);
for i = 1:NumNeurons
    Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
    sort( Lambda_Hat_Spikes );
    Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
    if Y > mean(Lambda_Spikes(Neuron_Type(i))+1)
        Neuron_Activated(i) = 1;
    end
end

save('Pt4.mat');


%% Functions

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random(:,:) < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end