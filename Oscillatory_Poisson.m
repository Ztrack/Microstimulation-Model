clearvars; clc;

clearvars; clc;
load('InitialConditionsFull.mat')
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)
%% FFT of sin wave
Lambda = 40;
dt = 1e-3;
t = 0:dt:1-dt;
x = sin(2*pi*Lambda*t);
figure; plot(x); title('40 Hz sin wave');
Y = fft(x,length(t));
Pyy = Y.*conj(Y)/length(t);
f = 1000/length(t)*(0:length(t)/2);
figure; plot(f,Pyy(1:(length(t)/2)+1))
title('PSD of a 40 Hz Sin Wave');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 100]);

%% FFT of non-oscillating neuron
Lambda = 40;
NumTrials = 1;
dt = 1e-3;
t = 0:dt:1-dt;

for i = 1:1000
    x = Simple_PoissonGen(Lambda, dt, NumTrials);
    Y = fft(x,length(t));
    Pyy1(i,:) = Y.*conj(Y)/length(t);
end
Pyy = mean(Pyy1);
f = 1000/length(t)*(0:length(t)/2);
figure; plot(f,Pyy(1:(length(t)/2)+1));
title('PSD of a Non-oscillatory Neuron');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 10]);

%% FFT of oscillating neuron

Lambda = 40;
NumTrials = 1;
dt = 1e-3;
t = 0:dt:1-dt;

for i = 1:1000
    x = Oscillatory_PoissonGen(Lambda, dt, NumTrials);
    Y = fft(x,length(t));
    Pyy1(i,:) = Y.*conj(Y)/length(t);
end
Pyy = mean(Pyy1);
f = 1000/length(t)*(0:length(t)/2);
figure; plot(f,Pyy(1:(length(t)/2)+1));
title('PSD of an oscillatory Neuron');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([0 100]);

%% FFT Mean & sum for all neurons

NumTrials = 1000;
dt = 1e-3;
t = 0:dt:1-dt;
f = 1000/length(t)*(0:length(t)/2);
Lambda = zeros(NumNeurons,1); Lambda_Hat = zeros(NumNeurons,1);
Lambda(neuron.excitatory) = 20;
Lambda(neuron.inhibitory) = 40;
Lambda_Hat(neuron.excitatory) = randi([20 60],size(neuron.excitatory));
Lambda_Hat(neuron.inhibitory) = randi([40 80],size(neuron.inhibitory));

for i = 1:NumNeurons
    if Neuron_Oscillatory_Type(i) == 3
        for n = 1 : NumTrials
            x = Oscillatory_PoissonGen_Spike_Train(Lambda(i), dt, 1);
            Y = fft(x,length(t));
            Pyy1(n,:) = Y.*conj(Y)/length(t);
        end
    else
        for n = 1 : NumTrials
            x = PoissonGen_Spike_Train(Lambda(i), dt, 1);
            Y = fft(x,length(t));
            Pyy1(n,:) = Y.*conj(Y)/length(t);
        end
    end
    Pyy2(i,:)= mean(Pyy1);
end
Pyy_mean_Lambda = mean(Pyy2);
Pyy_mean_Lambda2 = mean(Pyy2(Neuron_Inhibitory_Oscillatory,:));
Pyy_mean_Lambda3 = mean(Pyy2(Neuron_Inhibitory_Non_Oscillatory,:));
Pyy_mean_Lambda4 = mean(Pyy2(neuron.excitatory,:));

for i = 1:NumNeurons
    if Neuron_Oscillatory_Type(i) == 3
        for n = 1 : NumTrials
            x = Oscillatory_PoissonGen_Spike_Train(Lambda_Hat(i), dt, 1);
            Y = fft(x,length(t));
            Pyy1(n,:) = Y.*conj(Y)/length(t);
        end
    else
        for n = 1 : NumTrials
            x = PoissonGen_Spike_Train(Lambda_Hat(i), dt, 1);
            Y = fft(x,length(t));
            Pyy1(n,:) = Y.*conj(Y)/length(t);
        end
    end
    Pyy3(i,:)= mean(Pyy1);
end
Pyy_mean_Lambda_Hat = mean(Pyy3);
Pyy_mean_Lambda_Hat2 = mean(Pyy3(Neuron_Inhibitory_Oscillatory,:));
Pyy_mean_Lambda_Hat3 = mean(Pyy3(Neuron_Inhibitory_Non_Oscillatory,:));
Pyy_mean_Lambda_Hat4 = mean(Pyy3(neuron.excitatory,:));

figure; 
plot(f,Pyy_mean_Lambda(1:(length(t)/2)+1)); hold on;
plot(f,Pyy_mean_Lambda_Hat(1:(length(t)/2)+1),'--'); hold on;
title('PSD of Neural Population');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([20 100]); legend('Pre Stimulation','Post Stimulation');

figure; 
plot(f,Pyy_mean_Lambda2(1:(length(t)/2)+1)); hold on;
plot(f,Pyy_mean_Lambda_Hat2(1:(length(t)/2)+1),'--'); hold on;
title('PSD of OSC Inhibitory Neurons');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([20 100]); legend('Pre Stimulation','Post Stimulation');

figure; 
plot(f,Pyy_mean_Lambda3(1:(length(t)/2)+1)); hold on;
plot(f,Pyy_mean_Lambda_Hat3(1:(length(t)/2)+1),'--'); hold on;
title('PSD of Non-OSC Inhibitory Neurons');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([20 100]); legend('Pre Stimulation','Post Stimulation');

figure; 
plot(f,Pyy_mean_Lambda4(1:(length(t)/2)+1)); hold on;
plot(f,Pyy_mean_Lambda_Hat4(1:(length(t)/2)+1),'--'); hold on;
title('PSD of Excitatory Neurons');
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([20 80]); legend('Pre Stimulation','Post Stimulation');

%%
% Create a function showing how m and r0 changes. How microstimulation modulates the r0 and the m. 
% m is the change between the overall rate and the peak of the activity. m tells us the change in amplitude of the oscillation. 
% r0 tells us the baseline of the firing rate. 
% vary the microstimulation steps and compute r0 and m for each step

Lambda = 40;
NumTrials = 1;
dt = 1e-3;
t = 0:dt:1-dt;
Lambda_Vector = 40:0.1:80;

for i = 1:1000
    x = Oscillatory_PoissonGen_Spike_Train(Lambda, dt, NumTrials);
    Y = fft(x,length(t));
    Pyy1(i,:) = Y.*conj(Y)/length(t);
end
Pyy = mean(Pyy1);
Pyy(1:5) = Pyy(6); % Normalize first few Hz for easier analysis
f = 1000/length(t)*(0:length(t)/2);
r0 = median(Pyy); %rate 1 which corresponds to the amplitude of the FFT
m0 = max(Pyy); % gain which corresponds to the maximum of the frequency spectrum

for ii = 1:length(Lambda_Vector)
    for i = 1:1000
        x = Oscillatory_PoissonGen_Spike_Train(Lambda_Vector(ii), dt, NumTrials);
        Y = fft(x,length(t));
        Pyy1(i,:) = Y.*conj(Y)/length(t);
    end
    Pyy = mean(Pyy1);
    Pyy(1:5) = Pyy(6); % Normalize first few Hz for easier analysis
    f = 1000/length(t)*(0:length(t)/2);
    r(ii) = median(Pyy); %rate 1 which corresponds to the amplitude of the FFT
    m(ii) = max(Pyy); % gain which corresponds to the maximum of the frequency spectrum
end

figure; plot(Lambda_Vector,abs(m-m0)); title('Stimulus Modulation of Sinusoid Gain'); xlabel('Firing Rate (Hz)'); ylabel('Absolute Power Change');
figure; plot(Lambda_Vector,abs(r-r0)); title('Stimulus Modulation of Sinusoid Amplitude'); xlabel('Firing Rate (Hz)'); ylabel('Absolute Power Change');


%% Functions
function Spikes = PoissonGen_Spike_Train(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
end

function Spikes = Oscillatory_PoissonGen_Spike_Train(freqOscillation, dt, NumTrials)
dur = 1;
modIndex = 1;
NumSteps = 1/dt;
t = 1:1:NumSteps;
Spike_Probability = (freqOscillation * (1/dur) * (1 + modIndex .* cos(2*pi * freqOscillation * t/NumSteps))) * dt;
X_random = 1 + (0-1).*rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
end


function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end

function Trial_Spikes = Oscillatory_PoissonGen(freqOscillation, dt, NumTrials)
dur = 1;
modIndex = 1;
NumSteps = 1/dt;
t = 1:1:NumSteps;
Spike_Probability = (freqOscillation * (1/dur) * (1 + modIndex .* cos(2*pi * freqOscillation * t/NumSteps))) * dt;
X_random = 1 + (0-1).*rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end