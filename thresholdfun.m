clearvars; clc;
load('InitialConditionsFull.mat');
% Function to convert raw output to thresholds
% output = num electrodes x num steps x num neurons x num simulations
% output = 100 x 50 x 1000 x 2000 in typical settings

% Setup parameters
bpct = 05; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
NumTrials = 50; % Number of trials per poisson process
simulation = 2000; % Number of overall repeats for monte carle simulation
neuron.lambda = zeros(1,NumNeurons);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
inhibitoryfactor = 0.01; % at rate = 40hz (for inhibitory), there is a X% inhibition factor active. This increases linearly with lambda.
neuron.lambdamod(neuron.type == 1) = 40;
neuron.lambdamod(neuron.type == 2) = neuron.lambda(2) - neuron.lambda(1)*(neuron.lambda(1).*(inhibitoryfactor(1)/40));


%% Parsing large output data

% Must be done one at a time due to size of matrices
load('SEoutputc.mat');
for i = 1:size(output,1) % For every electrode
    for ii = 1:size(output,2)
        for iii = 1:size(output,3)
            outputbinary.c(i,ii,iii) = sum(output(i,ii,iii,:) > neuron.lambdamod(iii)+1) >= size(output,4).*.50;
        end
    end
end
clear('output');

load('SEoutputo.mat');
for i = 1:size(output,1) % For every electrode
    for ii = 1:size(output,2)
        for iii = 1:size(output,3)
            outputbinary.o(i,ii,iii) = sum(output(i,ii,iii,:) > neuron.lambdamod(iii)+1) >= size(output,4).*.50;
        end
    end
end

clear('output');
save('outputbinary.mat','outputbinary');

%% Constructing theshold per electrode
% We must find a proper threshold for each elecrode. This will be based
% on the first x neurons of the motion orientation we are interested in
% new output must be a low and a high value per stimulation strategy (4 total)

x = 5; % The current/opto step must activate this many motion tuned neurons to store the step

for i = 1:size(outputbinary.c,1)
    
        for ii = 1:size(outputbinary.c,2)
            numactivated(1,ii) = sum(outputbinary.c(i,ii,neuron.motion.number)) >= x;
            numactivated(2,ii) = sum(outputbinary.o(i,ii,neuron.nonmotion.number)) >= x;
        end
        
        threshold.c(i) = find(numactivated(1,:) == 1,1);
        threshold.o(i) = find(numactivated(2,:) == 1,1);
end

save('threshold.mat','threshold');

