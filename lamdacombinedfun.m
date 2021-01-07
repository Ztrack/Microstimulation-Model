function [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,lambdatype)

% Calculate firing rates frc,fro
[frc,fro,fro2] = fifun(neuron,Ie_Neurons,Il_Neurons,lambdatype);

% lambdatype:
% 1 = MS only
% 2 = opto only
% 3 = MS + Optogenetics (all excitatory - affects all cells)
% indiscriminately)
% 4 = MS + Optogenetics (all silencing - affects all cells indiscriminately
% 5 = MS only, subthreshold
% 6 = opto only, subthreshold
% 7. MS (subthreshold) + Opto (+ Subthreshold) + opto (- supathresh)

% Calculating Lambda hat based off microstimulation
if lambdatype == 1
    lambdahat = neuron.lambda + frc; % MS only, supathreshold
elseif lambdatype == 2
    lambdahat = neuron.lambda + fro; % Opto only excitation, supathreshold
elseif lambdatype == 3
    lambdahat = neuron.lambda + frc + fro; % MS + opto for all, subthreshold
elseif lambdatype == 4
    lambdahat = neuron.lambda + frc - fro; % MS - opto for all, supathreshold
elseif lambdatype == 5
    lambdahat = neuron.lambda + frc; % MS
elseif lambdatype == 6
    lambdahat = neuron.lambda + fro; % MS 
elseif lambdatype == 7
    lambdahat = neuron.lambda + frc + fro - fro2; % MS + opto for all, subthreshold
end

% Ensuring a firing rate limit is applied
limit = 300;
lambdahat(lambdahat>limit) = limit;
lambdahat(lambdahat<0) = 0;

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
end

% Applying EPSP Effect from pyramidal neurons onto pyramidal neurons (And
% motion to motion)
lambdahatpreconnections = lambdahat; % Lambdahat unmodified by synaptic connections

lambdahat(neuron.connections(:,2)) = lambdahat(neuron.connections(:,2)) + lambdahat(neuron.connections(:,1)).*neuron.connections(:,4);
lambdahat(neuron.connections((neuron.connections(:,3) == 1),1)) = lambdahat(neuron.connections((neuron.connections(:,3) == 1),1)) + lambdahatpreconnections(neuron.connections((neuron.connections(:,3) == 1),2)).*neuron.connections((neuron.connections(:,3) == 1),4);


% Ensuring a firing rate limit is applied
%lambdahat(lambdahat>limit) = limit;
lambdahat(lambdahat<0) = 0;

end