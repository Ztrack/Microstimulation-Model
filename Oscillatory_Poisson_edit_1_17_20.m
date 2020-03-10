close all; clear all; clc
freqOscillation = [35:45];

fs = 1000; % Hz
dur = .67; % in seconds

dt = 1/fs;
NumSteps = dur * fs;
NumTrials = 2000;

Trial_Spikes = nan(NumTrials,1);
Trial_Spikes_2 = Trial_Spikes;

t = 1:1:NumSteps;

modIndex = 1;

figure
hold on

fr = nan(NumTrials,1);
spikesAll = zeros(length(t), NumTrials);

for n = 1 : NumTrials
    
    f = randi([1 length(freqOscillation)],1,1);
    
    fr(n) = freqOscillation(f);
    osc_lambda = (freqOscillation(f) * (1/dur) * (1 + modIndex .* cos(2*pi * freqOscillation(f) * t/NumSteps))) * dt;
    
    signal = 1 + (0-1).*rand(NumSteps,1);

    spikes = (signal < (osc_lambda'));
    Trial_Spikes(n) = sum(spikes);
    spikesAll(:,n) = spikes;
    plot(t,spikes*n,'.b');

end


figure; plot(Trial_Spikes,'b');

% close all;
% for ff = 1 : length(freqOscillation)
%     figure; hold on
%     h = find(fr == freqOscillation(ff));
%     for n = 1 : length(h)
%         plot(t,spikesAll(:,h(n))*n,'.b');
%     end
%     title(num2str(freqOscillation(ff)))
%     pause
% end
% 
% 
