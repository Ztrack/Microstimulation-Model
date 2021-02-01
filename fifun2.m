function [frc,fro,fro2] = fifun(neuron,params,Ie_Neurons,Il_Neurons,lambdatype)

% y1 = Firing rate change due to current
% y2 = Firing rate change due to Optogenetics

% lambdatype:
% 1 = MS only, supathreshold
% 2 = opto only, supathreshold
% 3 = MS + Optogenetics (all excitatory - affects all cells)
% indiscriminately), both subthreshold
% 4 = MS + Optogenetics (all silencing - affects all cells
% indiscriminately, both supathreshold
% 5 = MS only, subthreshold
% 6 = opto only, subthreshold
% 7. MS + Opto (+ Subthreshold) + opto (- supathresh)

%% Current
frc = zeros(length(Ie_Neurons),1); % Initialize firing rate due to current

% Supathreshold case
frc(neuron.type == 1) = 0.40619.*(Ie_Neurons(neuron.type == 1)-0.22253);
frc(neuron.type == 2) = 0.140.*(Ie_Neurons(neuron.type == 2)-0.186);
frc(frc<0) = 0; % Fix negative

%Subthreshold Case
if lambdatype == 3 | lambdatype == 5 | lambdatype == 7 % Subthreshold case
        for i = 1:length(frc)
            if frc(i) == 0 % Did not reach threshold
                frc(i) = Ie_Neurons(i)/0.22253;
                frc(i) = Ie_Neurons(i)/0.186;
            end
        end
end

% Plot test:
% I0 = 0:1:500; ElectroNo = 45; 
% Ie_Neurons = neuron.io.soma(2,ElectrodeNo).*neuron.dirmult(1,ElectrodeNo).*I0 + neuron.io.axon(2,ElectrodeNo).*I0;
% frc = 0.140.*(Ie_Neurons-0.186); plot(Ie_Neurons,frc);
%% Opto
fro = zeros(length(Ie_Neurons),1);
fro2 = zeros(length(Ie_Neurons),1);
n = 0.82; % Hill coefficient
Imax = 26.6; % Maximum FR
k = 0.49; % half maximal light sensitivity
if lambdatype == 3 | lambdatype == 6 % Subthreshold case active
    fro(Il_Neurons < .10) = (Il_Neurons(Il_Neurons < .10)/.10);
    fro(Il_Neurons >= .10) = Imax .* ((Il_Neurons(Il_Neurons >= .10).^n)./((k.^n)+(Il_Neurons(Il_Neurons >= .10).^n)));
else % Supathreshold
    fro = Imax .* ((Il_Neurons.^n)./((k.^n)+(Il_Neurons.^n)));
end

if lambdatype == 7 % Subthreshold case active
    fro2(Il_Neurons < .10) = (Il_Neurons(Il_Neurons < .10)/.10);
    fro2(Il_Neurons >= .10) = Imax .* ((Il_Neurons(Il_Neurons >= .10).^n)./((k.^n)+(Il_Neurons(Il_Neurons >= .10).^n)));
end

% Plot test:
% I0 = 0:0.1:0.1; ElectroNo = 45; 
% Il_Neurons = neuron.oo.soma(2,ElectrodeNo).*I0; 
% fro = Imax .* ((Il_Neurons.^n)./((k.^n)+(Il_Neurons.^n))); plot(Il_Neurons,fro);
end