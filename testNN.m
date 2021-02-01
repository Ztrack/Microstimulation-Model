% Synaptic model integrating received spikes, and a probabilistic spiking model
% Source: https://www.jneurosci.org/content/32/29/9817
clearvars; clc;

dt = 1/1000;
numneurons = 900;
EIratio = 0.80;
C = 1; % constant determining the connectivity probability.
T_total = 1000; % Time in ms for experiment length, number of time bins
T_steps = 1; % Step length in ms

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
WEE = 1/((72.36-41.40)/0.65); % Excitatory-Excitatory synaptic weight
WEI = -1/((72.36-41.40)/3.48); % Excitatory - Inhibitory synaptic weight 
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
tI = 10; % ms time constant for synaptic decay
dIdt = zeros(numneurons,1); % Initialize Synaptic decay
dPdt = zeros(numneurons,1); % Initialize Synaptic decay
Ps = zeros(numneurons,1); % Initialize probability of spiking
raster = zeros(1000,length(S));
nbins = T_total/T_steps;

% Neuron type data taken from https://www.neuroelectro.org/
tP = zeros(numneurons,1); tP(neuron_type == 1) = 9.57; tP(neuron_type == 2) = 19.73; % Membrane Time constant (ms)- Time constant for the membrane to repolarize after a small current injection of fixed amplitude and duration
Pr = zeros(numneurons,1); Pr(neuron_type == 1) = 40/nbins; Pr(neuron_type == 2) = 20/nbins; % Background spiking probability (Hz) - AP discharge rate in the absence of current injection or a stimulus
%lambdahat = randi([0 10],numneurons,1)/bins; % Random current injection
lambdahat = (zeros(numneurons,1)+0)/nbins;
Ps = Pr + lambdahat;

for t = 1:T_steps:T_total
    % Each time step (dt) of 1 ms starts with each neuron, i, updating Ii with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.
    for i = 1:numneurons
        
        % Update Synaptic Weights
        I(i) = I(i) + (sum(neuron_weights(i,:)'.*S)/nbins); % if previous bin spikes, S=1 applying a synaptic effect from that neuron.
    
        % Probability updating
        dIdt(i) = (0 - I(i))/tP(i); % Spiking decay
        
    end
    
    % Determine Spiking
    S = Ps + I > rand(numneurons,1); % We determine whether the neuron spikes with the probability PS and update the spiking vector for the next time step
   
    I = I + dIdt; % Apply spiking decay
    Ps(S == 1) = Pr(S == 1) + lambdahat(S == 1); %If a neuron spikes, the probability is reset to the reset value
    I(S == 1) = 0; % Reset synaptic strength if spike occured
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

figure;
x = 1:2;
sumspikes = sum(raster,1);
y = [mean(sum(raster(:,neuron_type == 1),1)) mean(sum(raster(:,neuron_type == 2),1))];
ystd = [std(sum(raster(:,neuron_type == 1),1)) std(sum(raster(:,neuron_type == 2),1))];
y95 = 1.96*[ystd(1)/sqrt(sum(neuron_type==1)) ystd(2)/sqrt(sum(neuron_type==2))];
bar(x,y); hold on;
er = errorbar(x,y,ystd,ystd);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
title('Population Firing Rate');
ylabel('firing rate (Hz)');
xlabel('Inhibitory                Excitatory');
ylim([0 max(y)+10]);

