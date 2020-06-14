function output = Simple_PoissonGen3(lambdas, dt, NumTrials,bpct)
NumSteps = 1/dt;
Spike_Probability = lambdas.*dt;
Spikes = rand(NumNeurons,NumTrials,NumSteps);
for i = 1:length(lambdas)
    Spikes2 = (Spikes(i,:,:) < Spike_Probability(i));
end
SumSpikes = sum(Spikes, 3);

output = prctile(SumSpikes',bpct);
end