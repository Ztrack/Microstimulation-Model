function [frc,fro,fro2] = fifun(neuron,params,Ie_Neurons,Il_Neurons,lambdatype)

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

a1 = 31.2; % Inhibitory slope
a2 = 29.94; % Excitatory Slope
thresh1 = 0.162; % Inhibitory Threshold
thresh2 = 2.2582; % Excitatory Threshold

for i = 1:length(neuron.type)
    
    % General Case
    if neuron.type(i) == 1 & Ie_Neurons(i) > thresh1 % Inhibitory type
        frc(i) = a1.*sqrt(Ie_Neurons(i)-.162);
    elseif neuron.type(i) == 2 & Ie_Neurons(i) > thresh2 % Excitatory type
        frc(i) = a2.*sqrt(Ie_Neurons(i)-thresh2);
    
    % Subthreshold Cases       
    elseif neuron.type(i) == 1 & Ie_Neurons(i) < thresh1 & any(lambdatype == [3 5 7]) % Inhibitory type, less than threshold
        frc(i) = Ie_Neurons(i)/thresh1;
    elseif neuron.type(i) == 2 & Ie_Neurons(i) < thresh2 & any(lambdatype == [3 5 7])% Excitatory type, less than threshold
        frc(i) = Ie_Neurons(i)/thresh2;
    end
end

%% Luminous Intensity
% High-speed mapping of synaptic connectivity using photostimulation in Channelrhodopsin-2 transgenic mice
% H. Wang, et al. 2007

fro = zeros(length(Il_Neurons),1);
fro2 = zeros(length(Il_Neurons),1);

% Hill Equation Parameters
n = 0.82; % Hill coefficient
k = 0.49; % half maximal light sensitivity
Imax = 26.6; % Maximum FR
shift_x = -.05739; % Threshold value for inhibitory neuron
shift_y = 1.042;
shift_x2 = -0.8; % Threshold value for excitatory neuron
shift_y2 = 1;

for i = 1:length(Il_Neurons)
    
    % General Supathreshold Case
    if neuron.type(i) == 1 & Il_Neurons(i) > abs(shift_x) % Inhibitory type
        fro(i) = Imax .* (((Il_Neurons(i)+shift_x).^n)./((k.^n)+((Il_Neurons(i)+shift_x).^n)))*shift_y;
    elseif neuron.type(i) == 2 & Il_Neurons(i) > abs(shift_x2) % Excitatory type
        fro(i) = Imax .* (((Il_Neurons(i)+shift_x2).^n)./((k.^n)+((Il_Neurons(i)+shift_x2).^n)))*shift_y2;
    
    % Subthreshold Cases       
    elseif neuron.type(i) == 1 & Il_Neurons(i) < abs(shift_x) & any(lambdatype == [3 6 7]) % Inhibitory type, less than threshold
        fro(i) = Il_Neurons(i)/shift_x;
    elseif neuron.type(i) == 2 & Il_Neurons(i) < abs(shift_x2) & any(lambdatype == [3 6 7])% Excitatory type, less than threshold
        fro(i) = Il_Neurons(i)/shift_x2;
    end
    
end

%% Sanity check 1:
% Il_Neurons = 0.01:0.01:10;
% n = 0.82; % Hill coefficient
% k = 0.49; % half maximal light sensitivity
% Imax = 26.6; % Maximum FR
% 
% % Example Inhibitory neuron
% shift_x = -.05739;
% shift_y = 1.042;
% fro = Imax .* (((Il_Neurons+shift_x).^n)./((k.^n)+((Il_Neurons+shift_x).^n)))*shift_y;
% fro(Il_Neurons < abs(shift_x)) = nan;
% 
% % Example Excitatory Neuron
% shift_x = -0.8;
% shift_y = 1;
% fro2 = Imax .* (((Il_Neurons+shift_x).^n)./((k.^n)+((Il_Neurons+shift_x).^n)))*shift_y;
% fro2(Il_Neurons < abs(shift_x)) = nan;
% 
% figure; 
% plot(Il_Neurons,fro); hold on; plot(Il_Neurons,fro2); 
% legend('inhib','excit'); xlabel('Luminance mW/mm^2'); ylabel('Firing Rate (Hz)');
% 
% %% Sanity check 2:
% figure;
% x = 0:0.002:3;
% Imax = 26.6; % Maximum FR
% fro = Imax .* ((x.^n)./((k.^n)+(x.^n))); % Original hill eq and parameters
% plot(x,fro);
% hold on;
% 
% % Adjusted parameters for excitatory neuron
% n = n;
% k = k;
% fro = Imax .* (((x-.8).^n)./((k.^n)+((x-.8).^n)));
% fro(x<0.8) = 0;
% plot(x,fro);
% 
% % Adjusted Parameters for inhibitory neuron
% fro = Imax .* (((x-.05739).^n)./((k.^n)+((x-.05739).^n)))*1.042;
% fro(x<0.05739) = 0;
% plot(x,fro);
% legend('Original Hill Eq','Adjusted Excitatory','Adjusted Inhibitory'); xlabel('Luminance mW/mm2'); ylabel('Firing Rate');



end