clearvars; clc;
load('InitialConditionsFull.mat');
% Function to convert raw output to thresholds
% output = num electrodes x num steps x num neurons x num simulations
% output = 100 x 50 x 1000 x 2000 in typical settings

%% Constructing theshold per electrode
% We must find a proper threshold for each elecrode. This will be based
% on the first x neurons of the motion orientation we are interested in
% new output must be a low and a high value per stimulation strategy (4 total)

x = 50; % The current/opto step must activate this many motion tuned neurons to store the step
load('SEoutputc.mat');
for i = 1:size(output,1) % For every electrode
    numactivated = zeros(1,size(output,2));
    
    for ii = 1:size(output,2) % For every step
        numactivated(ii) = sum(output(i,ii,:)) > x;
    end
    
    threshold.c(i) = units(find(numactivated == 1,1));
end

load('SEoutputo.mat');
for i = 1:size(output,1) % For every electrode
    numactivated = zeros(1,size(output,2));
    
    for ii = 1:size(output,2) % For every step
        numactivated(ii) = sum(output(i,ii,:)) > x;
    end
    
    threshold.o(i) = units(find(numactivated == 1,1));
end

save('threshold.mat','threshold');

