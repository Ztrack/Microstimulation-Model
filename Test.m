Lambda = 28:.1:30;
Lambda_Base = 19.6;
NumTrials = 10000;
dt = 1/1000;
Numrepeats = 100;
lambda_needed = nan(Numrepeats,1);

for ii = 1:Numrepeats
for i = 1:length(Lambda)
    if isnan(lambda_needed(ii))
    Trial_Spikes = Simple_PoissonGen(Lambda(i), dt, NumTrials);
    Y = prctile(Trial_Spikes,05);
    if Y > Lambda_Base+1
        lambda_needed(ii,1) = Lambda(i);
        break
    end
    end
end
disp(ii)
end
optimized_lambda = mean(lambda_needed);


function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end