function [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,inhibitoryfactor,lambdatype)

% Calculate firing rates frc,fro
[frc,fro] = fifun(neuron,Ie_Neurons,Il_Neurons);

% lambdatype:
% 1 = MS only
% 2 = opto only
% 3 = MS + Optogenetics (all excitatory - affects all cells)
% indiscriminately)
% 4 = MS + Optogenetics (all silencing - affects all cells indiscriminately),
% 5 = MS + Optogenetics (only express excitatory opsin in inhibitory neurons)
% 6 = MS + Optogenetics (express excitatory opsin in inhibitory neurons & express inhibitory opsin in excitatory neurons)
% 7 = MS + Optogenetics (Express excitatory opsin in excitatory cells +
% express inhibitory opsin in all cells indiscriminately)

% Calculating Lambda hat based off microstimulation
if lambdatype == 1
    lambdahat = neuron.lambda + frc; % MS only
elseif lambdatype == 2
    lambdahat = neuron.lambda + fro; % Opto only excitation
elseif lambdatype == 3
    lambdahat = neuron.lambda + frc + fro; % MS + opto for all
elseif lambdatype == 4
    lambdahat = neuron.lambda + frc - fro; % MS - opto for all
elseif lambdatype == 5
    lambdahat = neuron.lambda + frc; % MS 
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + fro(neuron.inhibitory); % Optogenetics excitatory opsin for inhibitory
elseif lambdatype == 6
    lambdahat = neuron.lambda + frc; % MS
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + fro(neuron.inhibitory); % Optogenetics excitatory opsin for inhibitory
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) - fro(neuron.excitatory); % Optogenetics inhibitory opsin for excitatory
elseif lambdatype == 7
    lambdahat = neuron.lambda + frc; % MS
    lambdahat = lambdahat - fro; % Inhibitory opsin in all cells
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) + fro(neuron.excitatory); % Excitatory opsin in excitatory neurons
end

% Ensuring a firing rate limit is applied
limit = 300;
lambdahat(lambdahat>limit) = limit;

% Calculating Inhibition effect on each motif. Rate-based calculation
% Effect = summation of (new-old)*factor for each inhibitory in motif
Inhibitory_Effect = zeros(1,length(neuron.motif));
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(neuron.motif(i)) = Inhibitory_Effect(i) + lambdahat(neuron.inhibitory(i)).*(inhibitoryfactor/neuron.lambda(1));
end

% Applying Inhibition onto excitatory neurons
% New lambda = last calculation - inhibitory effect of that motif
for i = 1:length(neuron.excitatory)
    lambdahat(neuron.excitatory(i)) = lambdahat(neuron.excitatory(i)) - lambdahat(neuron.excitatory(i)).*Inhibitory_Effect(neuron.motif(i));
end

end