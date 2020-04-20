function [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Ir_Neurons,inhibitoryfactor,lambdatype)

% lambdatype:
% 1 = MS + Optogenetics (all silencing),
% 2 = MS + Optogenetics (Inhibitory neuron excitation, only)
% 3 = MS + Optogenetics (Inhibitory neuron excitation, excitatory silencing)

% Calculating Lambda hat based off microstimulation
if lambdatype == 1
    lambdahat = neuron.lambda + Ie_Neurons' - Ir_Neurons';
elseif lambdatype == 2
    lambdahat = neuron.lambda + Ie_Neurons';
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + Ir_Neurons(neuron.inhibitory)';
elseif lambdatype == 3
    lambdahat = neuron.lambda + Ie_Neurons';
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + Ir_Neurons(neuron.inhibitory)';
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) - Ir_Neurons(neuron.inhibitory)';
end

% Calculating Inhibition effect on each motif. Rate-based calculation
% Effect = summation of (new-old)*factor for each inhibitory in motif
Inhibitory_Effect = zeros(1,length(neuron.motif));
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(neuron.motif(i)) = Inhibitory_Effect(i) + (lambdahat(neuron.inhibitory(i)) - neuron.lambda(neuron.inhibitory(i))) .* inhibitoryfactor;
end

% Applying Inhibition onto excitatory neurons
% New lambda = last calculation - inhibitory effect of that motif
for i = 1:length(neuron.excitatory)
    lambdahat(neuron.excitatory(i)) = lambdahat(neuron.excitatory(i)) - Inhibitory_Effect(neuron.motif(i));
end

end