clearvars;
load('InitialConditions_Pt2.mat');

%% Finding RB for each neuron-Electrode Pair


Neuron_RB = NaN(1,NumNeurons);

h = 50; % Step size
I0 = h:h:20000;
numrepeats = 100;

parfor jj = 1:numrepeats
Ie_Axon_Neurons = zeros(1,NumNeurons); % Stores Extracellular Current values at neuron points, array form
Ie_Soma_Neurons = zeros(1,NumNeurons);
Neuron_RB1 = NaN(1,NumNeurons);
for ii = 1:length(I0)

Ie_Axon_Neurons = I0_Axon_Neurons.*I0(ii);
Ie_Soma_Neurons = I0_Soma_Neurons.*I0(ii);
Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons;

dt = 1/1000; % step/bin time length 1ms
RS_Lambda = 20; % Regular Spiking Neurons, Excitory

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection

for i = 1:NumNeurons
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))); % Lambda hat firing rate based off stimulation
end


NumTrials = 1000;
RS_Lambda_Spikes = mean(Simple_PoissonGen(20, dt, NumTrials)); % Excitatory Lambda Spikes

for i = 1:NumNeurons
    if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
        Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
        sort( Lambda_Hat_Spikes );
        Y = prctile(Lambda_Hat_Spikes,05); % Calculates bottom 5th percentile
        if Y > mean(RS_Lambda_Spikes)+1
            Neuron_RB1(i) = I0(ii);    
        end
    end
end
end
Neuron_RB(jj,:) = Neuron_RB1;
end

%% Plots & Data Analysis
load('Pt2.mat');
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)

% Options for plotting
options = struct;
options.handle     = figure;
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

Neuron_Excited_Per_Step1 = zeros(numrepeats,length(I0)); % Calculates the number of neurons excited at every step of current
for i = 1:length(I0)
    Neuron_Excited_Per_Step1(:,i) = sum(Neuron_RB <=I0(i),2)*100/length(Neuron_RB);
end

NeuronStimDist = NaN(1,NumNeurons);
for i = 1:NumNeurons
    NeuronStimDist(i) = sqrt((NeuronX(i)-ElectrodeX).^2 + (NeuronY(i) - ElectrodeY).^2);
end

Neuron_RB_Diff = zeros(1,NumNeurons/2);
for i = 1:NumNeurons/2
    Neuron_RB_Diff(:,i) = mean(Neuron_RB(:,i),1)-mean(Neuron_RB(:,Axon_Away(i)),1);
end

%plot_areaerrorbar(Neuron_Excited_Per_Step1,options); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); 
%title('Activated Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons'); xlim([0 200]); ylim([0 100]);

Axon_Toward_StimDist = NeuronStimDist(Axon_Toward);
Axon_Toward_Neuron_RB = mean(Neuron_RB(:,Axon_Toward));
Axon_Away_StimDist = NeuronStimDist(Axon_Away);
Axon_Away_Neuron_RB = Neuron_RB(1,Axon_Away);
figure; 
plot(Axon_Toward_StimDist,Axon_Toward_Neuron_RB,'.');
hold on
plot(Axon_Away_StimDist,Axon_Away_Neuron_RB,'.','color','red'); legend('Axon_Toward','Axon_Away');
hold on
y1 = (2.125e-07).*Axon_Toward_StimDist.^3+(0.0006138).*Axon_Toward_StimDist.^2+(0.2781).*Axon_Toward_StimDist+(-14.3);
y2 = (3.739e-07).*Axon_Away_StimDist.^3+(0.0003382).*Axon_Away_StimDist.^2+(0.8189).*Axon_Away_StimDist+(-57.67);
plot(Axon_Toward_StimDist,y1,'.','color','blue');
plot(Axon_Away_StimDist,y2,'.','color','red');
title('Neuron RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Stimulus AU'); legend('Toward Field','Away','Polynomial Fit','Polynomial Fit'); ylim([0 max(max(Neuron_RB))]);
ylim([0 8000]); xlim([0 2000]);

figure; plot(NeuronStimDist(1:NumNeurons/2),Neuron_RB_Diff,'.'); title('Neuron RB Difference vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Stimulus AU');
ylim([-4000 4000]);
%% Functions

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random(:,:) < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end

function plot_areaerrorbar(data, options)
    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'std';
    end
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = mean(data,1);
    data_std  = std(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    end
    
    % Plotting the result
    figure(options.handle);
    x_vector = [options.x_axis', fliplr(options.x_axis')];
    patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    
    if options.legendswitch == 1
    legend(options.legend);
    end
    hold off;
    
end