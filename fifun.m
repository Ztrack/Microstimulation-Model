function [frc,fro] = fifun(neuron,Ie_Neurons,Il_Neurons)

% y1 = Firing rate change due to current
% y2 = Firing rate change due to Optogenetics

%% Current Intensity
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766297/

frc = zeros(length(Ie_Neurons),1); % Initialize firing rate due to current
for i = 1:length(neuron.lambda)
    if neuron.type(i) == 1 & Ie_Neurons(i) > 0.162 % Inhibitory type
        frc(i) = 31.2.*sqrt(Ie_Neurons(i)-.162);
    elseif neuron.type(i) == 2 & Ie_Neurons(i) > 2.2582 % Excitatory type
        frc(i) = 29.94.*sqrt(Ie_Neurons(i)-2.2582);
    end
end

%% Luminous Intensity
% High-speed mapping of synaptic connectivity using photostimulation in Channelrhodopsin-2 transgenic mice
% H. Wang, et al. 2007

fro = zeros(length(Ie_Neurons),1);
n = 0.82; % Hill coefficient
Imax = 26.6; % Maximum FR
k = 0.49; % half maximal light sensitivity


fro = Imax .* ((Il_Neurons.^n)./((k.^n)+(Il_Neurons.^n)));

% Sanity check:
% x = 0:0.002:.1;
% fro = Imax .* ((x.^n)./((k.^n)+(x.^n)));
% plot(x,fro);
% hold on;
% Imax = 25; % Maximum FR
% fro = Imax .* ((x.^n)./((k.^n)+(x.^n)));
% plot(x,fro);
% legend('Imax 300','Imax 26'); xlabel('Luminance mW/mm2'); ylabel('Firing Rate');
end