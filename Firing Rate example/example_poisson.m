clear all; close all; clc
numTrials = 100;
FS = 1000;
deltaT = 1/FS;
epochLength = 1000;
t = 1:FS/1000:1000;
maxFire = [10];
numNeurons = 100; length(maxFire);
spikes = zeros(epochLength,numTrials,numNeurons);
load('RateFunction_1000ms.mat')

% maxFire = 20;

count = 0;
spikes_avg = zeros(epochLength,numNeurons);
for j = 1 : length(maxFire)
    
%     rateFunction_c = maxFire(j).*ones(epochLength,1);
    rateFunction_c = rateFunction .* maxFire(j);
    
    for n = 1 : numTrials
        
        signal = 1 + (0-1).*rand(epochLength,numNeurons);        
        for NN = 1 : numNeurons     
            if NN < round(numNeurons/2)
                spikes(:,n,NN) = (signal(:,NN) < rateFunction_c*deltaT);
            else
                spikes(:,n,NN) = (signal(:,NN) < rateFunction_c/2*deltaT);
            end
        end
        
        add_spikes = sum(randi([0 1],5,1));
        rr = randperm(size(spikes,1));
        if add_spikes ~= 0
            rr = rr(1:add_spikes);
            spikes(rr,n,:) = 1;
        else
        end

    end
    spikes_avg(:,j) = squeeze(mean(spikes(:,:,j),2,'omitnan'));
    
end

spikes = squeeze(spikes);

% plot(spikes_avg, 'o')

figure;
hold on
for ii = 1 : numTrials
    plot(t, spikes(:,ii) * (ii),'.k');
end


neuralPopulation_matrix =  zeros(100, 100);
P = randi([1 100], 100,2);
for ii = 1 : length(P)
   
    curNeuron = []; curFireNeuron = [];
    curNeuron = spikes(:,:,ii);
    curFireNeuron = mean(sum(curNeuron,1),'omitnan');
    neuralPopulation_matrix(P(ii,1),P(ii,2)) = curFireNeuron;
    
end

figure;
imagesc(neuralPopulation_matrix);