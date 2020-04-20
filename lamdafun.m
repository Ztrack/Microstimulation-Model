function [lambdahat] = lamdafun(neuron,Ie_Soma_Axon_Neurons,kk)

% Calculating Lambda hat based off microstimulation
lambdahat = neuron.lambda + Ie_Soma_Axon_Neurons';

% Calculating Inhibition effect on each motif. Rate-based calculation
% Effect = summation of (new-old)*factor for each inhibitory in motif
Inhibitory_Effect = zeros(1,length(neuron.motif));
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(neuron.motif(i)) = Inhibitory_Effect(i) + (lambdahat(neuron.inhibitory(i)) - neuron.lambda(neuron.inhibitory(i))) .* neuron.inhibitoryfactor(kk);
end

% Applying Inhibition onto excitatory neurons
% New lambda = last calculation - inhibitory effect of that motif
for i = 1:length(neuron.excitatory)
    lambdahat(neuron.excitatory(i)) = lambdahat(neuron.excitatory(i)) - Inhibitory_Effect(neuron.motif(i));
end

% Neuron connectivity between motion neurons - currently unused code

%             for i = 1:length(Neuron_Connected) % lambdahat Increase for motion tuned pairs
%                 lambdahat(Neuron_Connected(i,2)) = lambdahat(Neuron_Connected(i,2)) + neuron.lambda(2).*Axonal_Mult;
%                 if Neuron_Connected(i,3) == 1 % Bi-Directional connections
%                     lambdahat(Neuron_Connected(i,1)) = lambdahat(Neuron_Connected(i,1)) + neuron.lambda(2).*Axonal_Mult;
%                 end
%             end

end