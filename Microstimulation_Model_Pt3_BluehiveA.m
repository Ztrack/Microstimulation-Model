clearvars; clc;
load('InitialConditionsFull.mat')
Oscillatory_Behavior = 1; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location
Motion_Axon = 0; % if set =1, enabled, connections between relevant motion neurons is established
Axonal_Mult = .01;

Inhibitory_Factor = [0.01 0.03 0.1 0.125 0.15 0.2];

if Directional_Current_Modifier == 0
    Directional_Current_Mult(:,:) = 1;
end

if Motion_Axon == 0
    I0_Motion_Neurons = zeros(NumNeurons,length(electrode.x));
end

%% Extracellular Current Matrix /Loop Start
h = 50; % Step size
I0 = h:h:25000;
numrepeats = 100;
ElectrodeDist = sqrt((sx/2-electrode.x).^2 + (sy/2-electrode.y).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode


for kk = 1:length(Inhibitory_Factor)
    
    Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes Lambda+1 spike
    
    parfor jj = 1:numrepeats
        
        Neuron_RB1 = NaN(1,NumNeurons);
        
        for ii = 1:length(I0)
            Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
            Ie_Soma_Neurons = zeros(NumNeurons,1);
            
            for i = 1:length(ElectrodeNo) % Summation of current for every neuron component by every electrode & its corresponding current
                Ie_Axon_Neurons = Ie_Axon_Neurons + I0_Axon_Neurons(:,ElectrodeNo(i)).*I0(ii);
                Ie_Soma_Neurons = Ie_Soma_Neurons + I0_Soma_Neurons(:,ElectrodeNo(i)).*I0(ii).*Directional_Current_Mult(:,ElectrodeNo);
            end
            
            Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons
            
            % Multiplier Matrices
            dt = 1/1000; % step/bin time length 1ms
            t_total= 0:dt:T-dt ;
            NumSteps = 1/dt; % Number of steps/bins in time t_total
            FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
            RS_Lambda = 20; % Regular Spiking Neurons, Excitory
            
            Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection
            Lambda = zeros(NumNeurons,1);
            
            for i = 1:NumNeurons
                if neuron.type(i) == 1 % If it is an FS (Inhibitory)
                    Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
                    Lambda(i) = FS_Lambda; % Lambda regular, no microstimulation
                end
            end
            
            Inhibitory_Effect = zeros(NumMotifs,1);
            Inhibitory_Effect_Lambda = zeros(NumMotifs,1);
            for i = 1:NumMotifs
                Inhibitory_Effect(i) = Lambda_Hat(neuron.inhibitory(i)).*Inhibitory_Factor(kk); % Calculates the inhibitory effect from FS firing on RS neurons per motif
                Inhibitory_Effect_Lambda = FS_Lambda*Inhibitory_Factor(kk);
            end
            
            for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
                if neuron.type(i) == 2
                    Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(neuron.motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
                    Lambda(i) = RS_Lambda - Inhibitory_Effect_Lambda; % Lambda regular, no microstimulation
                else
                    % Inhib Neurons do not inhib themselves
                end
            end
            
%             for i = 1:length(Neuron_Connected) % Lambda_Hat Increase for motion tuned pairs
%                 Lambda_Hat(Neuron_Connected(i,2)) = Lambda_Hat(Neuron_Connected(i,2)) + RS_Lambda.*Axonal_Mult;
%                 if Neuron_Connected(i,3) == 1 % Bi-Directional connections
%                     Lambda_Hat(Neuron_Connected(i,1)) = Lambda_Hat(Neuron_Connected(i,1)) + RS_Lambda.*Axonal_Mult;
%                 end
%             end
            
            % Finding RB for each neuron
            Lambda_Spikes = [0 0];
            Lambda_Spikes(1) = FS_Lambda; % Inhibitory lambda spikes
            Lambda_Spikes(2) = RS_Lambda; % Excitatory Lambda Spikes
            NumTrials = 1000;
            
            for i = 1:NumNeurons
                if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
                    if intersect(Neuron_Inhibitory_Oscillatory,i) == i & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
                        Lambda_Hat_Spikes = Oscillatory_PoissonGen(Lambda_Hat(i), dt, NumTrials);
                    else % Otherwise use the simple function:
                        Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
                    end
                    
                    Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
                    if Y > mean(Lambda_Spikes(neuron.type(i))+1)
                        Neuron_RB1(i) = I0(ii);
                    end
                end
            end
        end
        Neuron_RB(jj,:) = Neuron_RB1;
    end
    
    %% Export
    if kk == 1
        Neuron_RBA1 = Neuron_RB;
        save('A1.mat','Neuron_RBA1');
    elseif kk ==2
        Neuron_RBA2 = Neuron_RB;
        save('A2.mat','Neuron_RBA2');
    elseif kk ==3
        Neuron_RBA3 = Neuron_RB;
        save('A3.mat','Neuron_RBA3');
    elseif kk ==4
        Neuron_RBA4 = Neuron_RB;
        save('A4.mat','Neuron_RBA4');
    elseif kk ==5
        Neuron_RBA5 = Neuron_RB;
        save('A5.mat','Neuron_RBA5');
    elseif kk ==6
        Neuron_RBA6 = Neuron_RB;
        save('A6.mat','Neuron_RBA6');
    end
end
%% Functions

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end

function Trial_Spikes = Oscillatory_PoissonGen(Firing_Rate, dt, NumTrials)
dur = 1;
modIndex = 1; % Higher modulation = higher variability in sine wave
NumSteps = 1/dt; % Number of time steps
t = 1:1:NumSteps;  % Time vector
Freq = randi([35 45],NumTrials,1); % Gamma Hz Oscillation
Spike_Probability = (Firing_Rate * (1/dur) * (1 + modIndex .* cos(2*pi * Freq * t/NumSteps))) * dt; % Probability of spike occuring in time step
X_random = 1 + (0-1).*rand(NumTrials,NumSteps); % Random Number that spike probability must overcome
Spikes = (X_random < Spike_Probability); % 1 = spike, 0 = no spike
Trial_Spikes = sum(Spikes, 2); % Sum of all spikes per trial
end