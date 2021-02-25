% Synaptic model integrating received spikes, and a probabilistic spiking model
% Source: https://www.jneurosci.org/content/32/29/9817
clearvars; clc;

dt = 1;
numneurons = 900;
EIratio = 0.80;
C = 1; % constant determining the connectivity probability.
T_total = 1000; % Time in ms for experiment length, number of time bins
dt = 1; % Step length in ms
rp = 1; % Refractory period, in ms

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
WIE = 0; % Inhibitory-Excitatory synaptic weight in mV
WEE = 0.65; % Excitatory-Excitatory synaptic weight in mV
WEI = -3.48; % Excitatory - Inhibitory synaptic weight in mV
WII = 0; % Inhibitory-Inhibitory synaptic weight in mV

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
Is = zeros(numneurons,1); % Initialize synaptic input from other neurons. 0 at time start.
S = zeros(numneurons,1); % Initialize a binary spiking vector of the neurons spiking in previous time step
raster = zeros(1000,length(S));
bins = T_total/dt;

% Neuron type data taken from https://www.neuroelectro.org/
tm = zeros(numneurons,1); tm(neuron_type == 1) = 9.57e-3; tm(neuron_type == 2) = 19.73e-3; % Membrane Time constant (ms)- Time constant for the membrane to repolarize after a small current injection of fixed amplitude and duration
vrest = zeros(numneurons,1); vrest(neuron_type == 1) = -67.52e-3; vrest(neuron_type == 2) = -72.36e-3; % Resting potential in mV
vthresh = zeros(numneurons,1); vthresh(neuron_type == 1) = -40.25e-3; vthresh(neuron_type == 2) = -41.40e-3;
v = vrest; % Initialize membrane potential for each neuron
R = zeros(numneurons,1); R(neuron_type == 1) = 3e7; R(neuron_type == 2) = 3e7;

% Current Injection, I(i) = I_R + I_C
I = zeros(numneurons,1); % Initialize current injection for each neuron
%lambdahat = randi([1 3],numneurons,1)/bins; % Random current injection

Irandom = 1e-9 .* random('Normal', 0, 1.5, [numneurons, bins]); % Input in nA
%Irandom = zeros(numneurons,bins);

for t = 1:dt:T_total
    
    % Synaptic Input in voltage
    for i = 1:numneurons
        % Each time step (t) of 1 ms starts with each neuron, i, updating Is with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.
        Is(i) = sum(neuron_weights(i,:)'.*S); % if previous bin spikes, S=1 applying a synaptic effect from that neuron.
    end
    
    % Membrane Potential
    for i = 1:numneurons
        % Each time step (t) of 1 ms starts with each neuron, i, updating Is with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.
        dv = (dt/tm(i)) .* (-v(i)+vrest(i) + (I(i)+Irandom(i,t)).*R(i)); % Change in membrane potential for this time step
    end
    v = v + dv; % Update neuron voltage
    
    % Determine Spiking
    S = v > vthresh; % We determine whether the neuron spikes
    
    
    % Reset
    v(S == 1) = vrest(S == 1); %If a neuron spikes, the voltage is reset to the resting potential
    raster(t,:) = S;
    
    %% 2
    V_reset = -0.080; % -80mV
    V_e = -0.075; % -75mV
    V_th = -0.040; % -40mV
    Rm = 10e6; % membrane resistance
    tau_m = 10e-3; % membrane time constant
    
    
    dt = 0.0002;
    T = 0:dt:1; % 1 second simulation
    
    Vm(1) = V_reset;
    Im = 5e-9;
    
    for t=1:length(T)-1,
        if Vm(t) > V_th,
            Vm(t+1) = V_reset;
        else,
            Vm(t+1) = Vm(t) + dt * ( -(Vm(t) - V_e) + Im * Rm) / tau_m;
        end;
    end;
    % on reset, add 1 to SR vector for that neuron
    plot(T,Vm,'b-');
    xlabel('Time(s)');
    ylabel('Voltage (V)');
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

