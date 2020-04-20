function Trial_Spikes = Simple_PoissonGen(lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end