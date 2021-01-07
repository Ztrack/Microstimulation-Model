% Synaptic model integrating received spikes, and a probabilistic spiking model
% Source: https://www.jneurosci.org/content/32/29/9817
clearvars; clc;

dt = 1/1000;
numneurons = 2500;
EIratio = 0.80;
C = 1; % constant determining the connectivity probability.

% Square area grid
% Given a simple square grid with 2500 neurons, the grid will be 50x50 long
% for simplicity

neuron_connected = zeros(numneurons,numneurons);
for i = 1:numneurons
    % Need to find the distance r from each neuron to every other neuron.
    Area = zeros(sqrt(numneurons),sqrt(numneurons));
    Area(i) = 1;
    Ed = bwdist(Area); % Ecludian distance
    r = Ed(:);
    
    % Calculating the connectivity
    P = C.*exp(-r); % Probability of synaptic connection P(r) = Ce^-r, where C is a constant determining the connectivity probability.
    P(i) = 0; % Chance of self-connecting is 0.
    neuron_connected(i,:) = P > rand(size(P)); % Binary vector containing connectivity information
end

% Neuron type definition
neuron_type = ones(numneurons,1)+1; neuron_type(numneurons*EIratio+1:end) = 1; % Initialize neuron population. 1 = inhibitory, 2 = excitatory
perm = randperm(length(neuron_type));
neuron_type = neuron_type(perm); % Randomize Excitatory and inhibitory neuron locations

% Synaptic Weights
% The weights (Wji) of these connections were dependent on the type (excitatory or inhibitory) of the presynaptic (i) and postsynaptic (j) neuron
WIE = 0; % Inhibitory-Excitatory synaptic weight
WEE = 0.1; % Excitatory-Excitatory synaptic weight
WEI = -WIE*10; % Excitatory - Inhibitory synaptic weight 
WII = 0; % Inhibitory-Inhibitory synaptic weight;

% Synaptic Weight vector
% Now we connect the synaptic weight to each neuron connectivity
neuron_weights = zeros(numneurons,numneurons);
for i = 1:numneurons
    x = find(neuron_connected(i,:) == 1); % Connected neurons
    for ii = 1:length(x)
        if neuron_type(i) == 1 & neuron_type(x(ii)) == 1 % If postsynaptic neuron is inhibitory and pre-synaptic is inhibitory, apply WII
            neuron_weights(i,x(ii)) = WII;
        elseif neuron_type(i) == 1 & neuron_type(x(ii)) == 2 % Inhibitory - Excitatory
            neuron_weights(i,x(ii)) = WIE;
        elseif neuron_type(i) == 2 & neuron_type(x(ii)) == 1 % Excitatory - Inhibitory
            neuron_weights(i,x(ii)) = WEI;  
        elseif neuron_type(i) == 2 & neuron_type(x(ii)) == 2 % Excitatory - Excitatory
            neuron_weights(i,x(ii)) = WEE;
        end
    end
end

% Spike Model
I = zeros(numneurons,1); % Initialize synaptic input from other neurons. 0 at time start.
S = zeros(numneurons,1); % Initialize a binary spiking vector of the neurons spiking in previous time step
tI = 9; % ms time constant
dIdt = zeros(numneurons,1); % Initialize Synaptic decay
dPdt = zeros(numneurons,1); % Initialize Synaptic decay
Ps = zeros(numneurons,1); % Initialize probability of spiking
% Neuron type data taken from https://www.neuroelectro.org/
tP = zeros(numneurons,1); tP(neuron_type == 1) = 9.57; tP(neuron_type == 2) = 19.73; % initialize probability of spiking decay in Hz
Pr = zeros(numneurons,1); Pr(neuron_type == 1) = 5.32; Pr(neuron_type == 2) = 0.67;% probability of spiking reset in Hz
raster = zeros(1000,length(S));

Ps = Pr;
for t = 1:1:1000
    % Each time step (dt) of 1 ms starts with each neuron, i, updating Ii with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.
    for i = 1:numneurons
        I(i) = I(i) + sum(neuron_weights(i,:)'.*S); % if previous bin spikes, S=1 applying a synaptic effect from that neuron.
        dIdt(i) = (0 - I(i))/tI; % Synaptic decay with tI time factor
    end
    I = I + dIdt; % Apply synaptic decay
    
    % Probability-based spiking
    for i = 1:numneurons
        Ps(i) = Ps(i) + I(i); % update spiking probability with synaptic weight
        dPdt(i) = (Pr(i) - Ps(i))/tP(i); % Spiking decay
    end
    Ps = Ps + dPdt; % Apply spiking decay
    
    S = Ps > rand(numneurons,1)*100; % We determine whether the neuron spikes with the probability PS and update the spiking vector for the next time step
    Ps(S == 1) = 0; %If a neuron spikes, the probability is reset to the reset value
    I(S == 1) = 0; % Reset synaptic strength
    
    raster(t,:) = S;
end

figure;
plot(1:1000,sum(raster,2));
title('Population Activity');
ylabel('Number of Spikes');
xlabel('Time (ms)');

figure;
hist(sum(raster,1),25);
title('Neuron Activity Across Population');
ylabel('Number of Neurons');
xlabel('Firing Rate (Hz)');

figure;
for i = 1:numneurons
    ylim([1 200]);
    xlim([0 1000]);
    plot(1:1000,raster(:,i).*i,'.');
    hold on;
end
title('Raster Data');
ylabel('Neuron Number');
xlabel('Time (ms)');
