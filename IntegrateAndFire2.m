function out_FR = IntegrateAndFire2(neuron,params,Im)

% Spike Model
%Im = zeros(params.numneurons,1)+0.2e-9; Input Current
Vs = zeros(params.numneurons,1); % Initialize synaptic input from other neurons as a voltage. 0 at time start.
S = zeros(1,params.numneurons); % Initialize a binary spiking vector of the neurons spiking in previous time step
out_FR = zeros(params.numneurons,1); % Initialize output rate
dt = 0.0002; % Time step length
T = 0:dt:1; % 1 second simulation
%raster = zeros(1000,length(S));

% Neuron type data taken from https://www.neuroelectro.org/
In = zeros(params.numneurons,1); In(neuron.type == 1) = 0.1991e-9; In(neuron.type == 2) = 0.2165e-9; % Noisy current input
V_e = zeros(params.numneurons,1); V_e(neuron.type == 1) = -67.52e-3; V_e(neuron.type == 2) = -72.36e-3; % Resting Membrane Potential
V_th = zeros(params.numneurons,1); V_th(neuron.type == 1) = -40.25e-3; V_th(neuron.type == 2) = -41.40e-3; % Threshold Membrane Potential
V_reset = zeros(params.numneurons,1); V_reset(neuron.type == 1) = -80e-3; V_reset(neuron.type == 2) = -80e-3; % Membrane Potential Reset after spike
Rm = zeros(params.numneurons,1); Rm(neuron.type == 1) = 152.7e6; Rm(neuron.type == 2) = 160.47e6;
tau_m = zeros(params.numneurons,1); tau_m(neuron.type == 1) = 9.57e-3; tau_m(neuron.type == 2) = 19.73e-3; % Membrane Time constant (ms)- Time constant for the membrane to repolarize after a small current injection of fixed amplitude and duration

Vm = V_reset; % Initiate Membrane potential

for t=1:length(T)-1
    
    % Each time step (dt) of 1 ms starts with each neuron, i, updating Ii with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.
    % Update Synaptic Weights
    Vs = sum(neuron.weight.matrix.*S,2); % if previous bin spikes, S=1 applying a synaptic effect from that neuron.
    
    % Determine Spiking
    S = Vm > V_th; % We determine whether the neuron spikes by crossing threshold
    Vm(S == 1) = V_reset(S == 1); %If a neuron spikes, the membrane potential Vs reset to the reset value
    
    % Update membrane potential
    Vm = Vm + Vs + dt .* ( -(Vm - V_e) + In + Im .* Rm) ./ tau_m;
    
    % Output
    %raster(t,:) = S; % For debugging
    out_FR = out_FR + S;
end

end