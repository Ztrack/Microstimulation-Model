function output = nonhomoPoissonGen(lambda, dt, numTrials,simulation,bpct,rateFunction)
% New non-homogenous poisson funct, a no homo poisson if you will
% combining MGR's "non_homogenous_poisson().m",
% with "Simple_PoissonGen.m" from HB's Microstimulation_Model_Pt3_BluehiveA.m
% match input/output of Simple_PoissonGen.m
%
% SAMPLE INPUT
% (28,0.001,100)

%PREVIOUS CODE Simple_PoissonGen.m
%----------------------------
% %function Trial_Spikes = Simple_PoissonGen(lambda, dt, NumTrials)
% %   NumSteps = 1/dt;
% %   Spike_Probability = lambda.*dt;
% %   X_random = rand(NumTrials,NumSteps);
% %   Spikes = (X_random < Spike_Probability);
% %   Trial_Spikes = sum(Spikes, 2);
% %end
%
% NOTES
% 20200501 mduhain
%  - function created
% 20200610 mduhain
%  - binning added to output (10 total, 100ms each)

%% Load Vars
rateFunction = (rateFunction .* 0.9) + 3; %add baseline firing rate of 3 sp/s

%% MAIN
output = zeros(simulation,10);

for i = 1:simulation
    spikes = zeros(length(rateFunction),numTrials); %empty mat
    rasterPlot = spikes; %local copy
    scalar = lambda/14.34; %relates lambdaHat value, with this non-homogenous poisson's rate (avg 14.34)
    signal = 1 + (0-1).*rand(length(rateFunction),numTrials); %noise
    % MGR explained the rationale, but why when "1+(0-1).*noise" == "noise"?
    spikes = (signal < rateFunction'*dt*scalar);% 0,1 raster for spikes w/ noise, based on rateFunct
    %collapse collumns, so each row(trial) is the sum of all spikes present
    %Trial_Spikes = sum(spikes,1)';
    
    Alt_mat = reshape(spikes',[numTrials,100,10]);
    % each row represents data from a single neuron
    % each column represents boolean firing at each ms
    % each z-stack, i.e. 1st bin, all neurons x all trials (first 100ms)
    
    Trial_Spikes = squeeze(sum(Alt_mat,2));
    % rows represent trials
    % collums represent spikes per bin (100ms bin size);
    
    output(i,:) = prctile(Trial_Spikes,bpct); % Outputs bottom 5th percentile of spikes
    
end


end