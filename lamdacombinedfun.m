function [lambdahat,lambdamod] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,lambdatype)

% Calculate firing rates frc,fro
[frc,fro] = fifun(neuron,Ie_Neurons,Il_Neurons,lambdatype);

% lambdatype:
% 1 = MS only
% 2 = opto only
% 3 = MS + Optogenetics (all excitatory - affects all cells)
% indiscriminately)
% 4 = MS + Optogenetics (all silencing - affects all cells indiscriminately
% 5 = MS+Opto(+all 50%) - opto(-all not bounded)
% 6 = MS + Optogenetics (only express excitatory opsin in inhibitory neurons)
% 7 = MS + Optogenetics (express excitatory opsin in inhibitory neurons & express inhibitory opsin in excitatory neurons)
% 8 = MS + Optogenetics (Express excitatory opsin in excitatory cells +
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
    lambdahat = neuron.lambda + frc + fro - fro;
elseif lambdatype == 6
    lambdahat = neuron.lambda + frc; % MS 
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + fro(neuron.inhibitory); % Optogenetics excitatory opsin for inhibitory
elseif lambdatype == 7
    lambdahat = neuron.lambda + frc; % MS
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + fro(neuron.inhibitory); % Optogenetics excitatory opsin for inhibitory
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) - fro(neuron.excitatory); % Optogenetics inhibitory opsin for excitatory
elseif lambdatype == 8
    lambdahat = neuron.lambda + frc; % MS
    lambdahat = lambdahat - fro; % Inhibitory opsin in all cells
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) + fro(neuron.excitatory); % Excitatory opsin in excitatory neurons
    
end

% Ensuring a firing rate limit is applied
limit = 300;
lambdahat(lambdahat>limit) = limit;
lambdahat(lambdahat<0) = 0;

% Neuron Lambda Modified
lambdamod = neuron.lambda; % This is the baseline firing rate after EPSP/IPSP effects are calculated

% Calculating Inhibition effect on each motif. Rate-based calculation
% Effect = summation of (new-old)*factor for each inhibitory in motif
Inhibitory_Effect = zeros(1,max(neuron.motif));
for i = 1:length(neuron.inhibitory) % For every inhibitory neuron, an inhibitory effect onto the motif is created
    Inhibitory_Effect(neuron.motif(neuron.inhibitory(i))) = Inhibitory_Effect(neuron.motif(neuron.inhibitory(i))) + lambdahat(neuron.inhibitory(i)).*neuron.inhibitoryfactor(i);
end

% Applying Inhibition onto excitatory neurons
% New lambda = last calculation - inhibitory effect of that motif
for i = 1:length(neuron.excitatory)
    lambdahat(neuron.excitatory(i)) = lambdahat(neuron.excitatory(i)) - Inhibitory_Effect(neuron.motif(neuron.excitatory(i)));
    lambdamod(neuron.excitatory(i)) = lambdamod(neuron.excitatory(i)) - Inhibitory_Effect(neuron.motif(neuron.excitatory(i)));
end

% Applying EPSP Effect from pyramidal neurons onto pyramidal neurons (And
% motion to motion)
for i = 1:length(neuron.connections)
    lambdahat(neuron.connections(i,2)) = lambdahat(neuron.connections(i,2)) + lambdahat(neuron.connections(i,1))*neuron.connections(i,4);
    lambdamod(neuron.connections(i,2)) = lambdamod(neuron.connections(i,2)) + lambdamod(neuron.connections(i,1))*neuron.connections(i,4);
    if neuron.connections(i,3) == 1 % If bi-directional
        lambdahat(neuron.connections(i,1)) = lambdahat(neuron.connections(i,1)) + lambdahat(neuron.connections(i,2))*neuron.connections(i,4);
        lambdamod(neuron.connections(i,1)) = lambdamod(neuron.connections(i,1)) + lambdamod(neuron.connections(i,2))*neuron.connections(i,4);
    end
end


end