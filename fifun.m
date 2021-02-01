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


%% Current Intensity
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766297/

frc = zeros(length(Ie_Neurons),1); % Initialize firing rate due to current

% for i = 1:length(neuron.lambda)
%     if neuron.type(i) == 1 & any(lambdatype == [1 2 4]) & Ie_Neurons(i) > 0.162 % Inhibitory type
%         frc(i) = 31.2.*sqrt(Ie_Neurons(i)-.162);
%     elseif neuron.type(i) == 2 & any(lambdatype == [1 2 4]) & Ie_Neurons(i) > 2.2582 % Excitatory type
%         frc(i) = 29.94.*sqrt(Ie_Neurons(i)-2.2582);
%     
%     % Case 3 - subthreshold
%     elseif neuron.type(i) == 1 & lambdatype == 3 & Ie_Neurons(i) > 0.162 % Inhibitory type
%         frc(i) = 31.2.*sqrt(Ie_Neurons(i)-.162);
%     elseif neuron.type(i) == 2 & lambdatype == 3 & Ie_Neurons(i) > 2.2582 % Excitatory type
%         frc(i) = 29.94.*sqrt(Ie_Neurons(i)-2.2582);
%     elseif neuron.type(i) == 1 & lambdatype == 3 & Ie_Neurons(i) < 0.162 % Inhibitory type
%         frc(i) = (Ie_Neurons(i)/0.162);
%     elseif neuron.type(i) == 2 & lambdatype == 3 & Ie_Neurons(i) < 2.2582 % Excitatory type
%         frc(i) = (Ie_Neurons(i)/2.2582);
%     end
% end

% for i = 1:length(neuron.lambda)
%         if any(lambdatype == [3 5 7]) & Ie_Neurons(i) < 2.2582 % Subthreshold case, EPSP is active
%             frc(i) = (Ie_Neurons(i)/2.2582);
%         elseif Ie_Neurons(i) > 2.2582 % Supathreshold case, no EPSP
%             frc(i) = 29.94.*sqrt(Ie_Neurons(i)-2.2582);
%         end
% end

if lambdatype == 3 | lambdatype == 5 | lambdatype == 7 % Subthreshold case
    frc(Ie_Neurons < 2.2582) = Ie_Neurons(Ie_Neurons < 2.2582)/2.2582;
    frc(Ie_Neurons >= 2.2582) = 29.94.*sqrt(Ie_Neurons(Ie_Neurons >= 2.2582)-2.2582);
else % Supathreshold case
    frc = 29.94.*sqrt(Ie_Neurons-2.2582);
end
r = isreal(frc);

%% Luminous Intensity
% High-speed mapping of synaptic connectivity using photostimulation in Channelrhodopsin-2 transgenic mice
% H. Wang, et al. 2007

fro = zeros(length(Ie_Neurons),1);
fro2 = zeros(length(Ie_Neurons),1);
n = 0.82; % Hill coefficient
Imax = 26.6; % Maximum FR
k = 0.49; % half maximal light sensitivity
% fro = Imax .* ((Il_Neurons.^n)./((k.^n)+(Il_Neurons.^n)));

% Sanity check:
% x = 0:0.002:.1;
% fro = Imax .* ((x.^n)./((k.^n)+(x.^n)));
% plot(x,fro);
% hold on;
% Imax = 26.6; % Maximum FR
% fro = Imax .* ((x.^n)./((k.^n)+(x.^n)));
% plot(x,fro);
% legend('Imax 300','Imax 26'); xlabel('Luminance mW/mm2'); ylabel('Firing Rate');

% for i = 1:length(neuron.lambda)
%         if (any(lambdatype == [3 6])) & (Il_Neurons(i) < .10) % Subthreshold case - EPSP is active
%             fro(i) = (Il_Neurons(i)/.10);
%         elseif Il_Neurons(i) > .10 % Supathreshold
%             fro(i) = Imax .* ((Il_Neurons(i).^n)./((k.^n)+(Il_Neurons(i).^n)));
%         else
%             fro(i) = 0; % Did not reach threshold in the non-subthreshold case
%         end
% end
% 
% if lambdatype == 7
%     fro2 = zeros(length(Ie_Neurons),1); % Initiate second opsin
%     for i = 1:length(neuron.lambda)
%         if (Il_Neurons(i) < .10) % Subthreshold case active
%             fro2(i) = (Il_Neurons(i)/.10);
%         elseif Il_Neurons(i) > .10 % Supathreshold
%             fro2(i) = Imax .* ((Il_Neurons(i).^n)./((k.^n)+(Il_Neurons(i).^n)));
%         else
%             fro2(i) = 0; % Did not reach threshold in the non-subthreshold case
%         end
%     end
% end


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

end