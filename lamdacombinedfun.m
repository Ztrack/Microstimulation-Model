function [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Ir_Neurons,inhibitoryfactor,lambdatype)

% lambdatype:
% 1 = MS only
% 2 = opto only
% 3 = MS + Optogenetics (all excitatory - affects all cells
% indiscriminately)
% 4 = MS + Optogenetics (all silencing - affects all cells indiscriminately),
% 5 = MS + Optogenetics (only express excitatory opsin in inhibitory neurons)
% 6 = MS + Optogenetics (express excitatory opsin in inhibitory neurons & express inhibitory opsin in excitatory neurons)
% 7 = MS + Optogenetics (Express excitatory opsin in excitatory cells +
% express inhibitory opsin in all cells indiscriminately)

% Calculating Lambda hat based off microstimulation
if lambdatype == 1
    lambdahat = neuron.lambda + neuron.lambda.*Ie_Neurons'; % MS only
elseif lambdatype == 2
    lambdahat = neuron.lambda + neuron.lambda.*Ir_Neurons'; % Opto only excitation
elseif lambdatype == 3
    lambdahat = neuron.lambda + neuron.lambda.*(Ie_Neurons' + Ir_Neurons'); % MS + opto for all
    elseif lambdatype == 4
    lambdahat = neuron.lambda + neuron.lambda.*(Ie_Neurons' - Ir_Neurons'); % MS - opto for all
elseif lambdatype == 5
    lambdahat = neuron.lambda + neuron.lambda.*Ie_Neurons'; % MS 
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + Ir_Neurons(neuron.inhibitory)'; % Optogenetics excitatory opsin for inhibitory
elseif lambdatype == 6
    lambdahat = neuron.lambda + neuron.lambda.*Ie_Neurons';
    lambdahat(neuron.inhibitory) = lambdahat(neuron.inhibitory) + Ir_Neurons(neuron.inhibitory)'; % Optogenetics excitatory opsin for inhibitory
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) - Ir_Neurons(neuron.inhibitory)'; % Optogenetics inhibitory opsin for inhibitory
elseif lambdatype == 7
    lambdahat = neuron.lambda + neuron.lambda.*Ie_Neurons';
    lambdahat = lambdahat - Ir_Neurons'; % Inhibitory opsin in all cells
    lambdahat(neuron.excitatory) = lambdahat(neuron.excitatory) + Ir_Neurons(neuron.excitatory)'; % Excitatory opsin in excitatory neurons
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