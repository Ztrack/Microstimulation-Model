function output = Simple_PoissonGen2(lambda, dt, NumTrials,simulations,bpct)
NumSteps = 1/dt;
Spike_Probability = lambda.*dt;
Spikes = rand(simulations,NumTrials,NumSteps);
Spikes = (Spikes < Spike_Probability);
SumSpikes = sum(Spikes, 3);

output = prctile(SumSpikes',bpct);
end