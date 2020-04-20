function Trial_Spikes = Oscillatory_PoissonGen(lambda, dt, NumTrials)
dur = 1;
modIndex = 1; % Higher modulation = higher variability in sine wave
NumSteps = 1/dt; % Number of time steps
t = 1:1:NumSteps;  % Time vector
Freq = randi([35 45],NumTrials,1); % Gamma Hz Oscillation
Spike_Probability = (lambda * (1/dur) * (1 + modIndex .* cos(2*pi * Freq * t/NumSteps))) * dt; % Probability of spike occuring in time step
X_random = 1 + (0-1).*rand(NumTrials,NumSteps); % Random Number that spike probability must overcome
Spikes = (X_random < Spike_Probability); % 1 = spike, 0 = no spike
Trial_Spikes = sum(Spikes, 2); % Sum of all spikes per trial
end