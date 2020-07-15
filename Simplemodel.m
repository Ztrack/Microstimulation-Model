clearvars; clc;

%Run with 1 neurons on a map and 1 electrode. Current vs FR. Use linear relationship of FR = .5x;
%Must run it so we find 50% of current and opto. Use together to get 1 spike. Show that any less results in no spikes.

%% Initial Condition Generation

sx = 100;
sy = 100;
map = zeros(sx,sy);
electrode.radii = 10;
electrode.x = size(map,1)/2; electrode.y = size(map,1)/2;

NumNeurons = 1;
neuron.type(1) = 1;
neuron.radii = 5;
neuron.x = 25; neuron.y = 25;

neuron.direction = randi([1,4],1,1); % 1 = x value down (left),2 = Y value down (up), 3 = X value up (right), 4 = Y value up (down)
re = [0,0,0,0]; re(neuron.direction) = 1; % Random element selection to choose orientation of axon
rl = randi([50+neuron.radii,480+neuron.radii],1); % Random length of axon + neuron.radii
xc = [neuron.x+(re(3)*neuron.radii)-(re(1)*neuron.radii) neuron.x+(re(3)*rl)-(re(1)*rl)]; % x coordinates of neuron, axon point
yc = [neuron.y+(re(4)*neuron.radii)-(re(2)*neuron.radii) neuron.y+(re(4)*rl)-(re(2)*rl)]; % y coordinates of neuron, axon point
if xc(2) <= 0 % Simple filter Makes sure the axon end does not go out of bounds in x
    xc(2) = 1;
elseif xc(2) >= sx
    xc(2) = sx;
end
Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
lindex = sub2ind([sx sy], rIndex,cIndex); % Linear indices
neuron.indices(1).axon = lindex(:);

for i = 1:NumNeurons
    singlepop = zeros(sx,sy);
    singlepop(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii, neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = 1;
    neuron.indices(i).soma = find(singlepop == 1);
end

Stim_Distance_Map = zeros(length(electrode.x),sx,sy);
for i = 1:length(electrode.x)
    Stim_Loc =  zeros(sx, sy); % Location of stimulus
    Stim_Loc(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
    Stim_Distance_Map(i,:,:) = Stim_Distance;
end
load('lightspread.mat');
neuron.io.soma = zeros(NumNeurons,length(electrode.x)); % Stores 1/r^2 area component of Somas
neuron.io.axon = zeros(NumNeurons,length(electrode.x)); % Stores 1/r^2 area component of Axons
neuron.oo.soma = zeros(NumNeurons,length(electrode.x)); % Stores light stimulus area summation component of Somas
lightspread.calc = [];

for i = 1:NumNeurons
    
    for j = 1:length(electrode.x)
        
        neuron.io.axon(i,j) = sum(sum(1./Stim_Distance_Map(j,neuron.indices(i).axon).^2)); % Summation 1/r^2 area component of axon
        neuron.io.soma(i,j) = sum(sum(1./Stim_Distance_Map(j,neuron.indices(i).soma).^2)); % Summation 1/r^2 area component of soma % Summation 1/r^2 area component of soma
        
        lightspread.calc = lightspread.averaged.a.*exp(lightspread.averaged.b.*Stim_Distance_Map(j,neuron.indices(i).soma)); % Light Intensity Calculation
        % lightspread.calc(lightspread.calc < 0) = 0; % Should not happen, debugging
        neuron.oo.soma(i,j) = sum(sum(lightspread.calc));
        
    end
    
end

population.neuron.map = zeros(sx,sy); population.axon.map = population.neuron.map; population.electrode.map = population.axon.map;% population of all elements
population.neuron.map(neuron.indices(1).soma) = 1;
population.neuron.map(neuron.indices(1).soma) = 1;
population.axon.map(neuron.indices(1).axon) = 1;
population.electrode.map(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1;

%% Lambdahat Calculation

bpct = 05; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
NumTrials = 100;
numrepeats = 100; % Number of overall repeats
ElectrodeNo = 1;
h = 100; % Steps
inhibitoryfactor = 0.01;

for kkk = 1:2
    lambdatype = kkk;
    
    if lambdatype == 1
        unitsmax = 100;
    else
        unitsmax = .05;
    end
    I0 = linspace(0,unitsmax,h);
    
    Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
    
    for jj = 1:numrepeats
        
        Neuron_RB1 = NaN(1,NumNeurons);
        
        for ii = 1:length(I0)
            
            Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*I0(ii) + neuron.io.axon(:,ElectrodeNo).*I0(ii); % Summation of current directly from stimulus + backpropogated up by axons. AU Current
            Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*I0(ii); % Summation of current directly from stimulus. AU irridance
            
            % Calculate neuron.lambda Change for poisson process
            frc = Ie_Neurons.*.5;
            fro = Il_Neurons.*.5;
            
            if lambdatype == 1
                lambdahat = neuron.lambda + frc; % MS only
            elseif lambdatype == 2
                lambdahat = neuron.lambda + fro; % Opto only excitation
            end
            limit = 300;
            lambdahat(lambdahat>limit) = limit;
            % Finding RB for each neuron
            dt = 1/1000;
            for i = 1:NumNeurons
                if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
                    Lambda_Hat_Spikes = Simple_PoissonGen(lambdahat(i), dt, NumTrials);
                    
                    Y = prctile(Lambda_Hat_Spikes,bpct); % Calculates bottom xth percentile
                    if Y > neuron.lambda(i)+1
                        Neuron_RB1(i,:) = I0(ii);
                    end
                end
            end
        end
        Neuron_RB(jj,:) = Neuron_RB1;
    end
    if lambdatype == 1
        solrb.e1 = Neuron_RB;
        solrb.e.I0 = I0;
    else
        solrb.o1 = Neuron_RB;
        solrb.o.I0 = I0;
    end
    
end

%% Analysis

stepsol.current.all = zeros(numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
stepsol.current.excitatory = stepsol.current.all; stepsol.current.inhibitory = stepsol.current.all; stepsol.current.motion = stepsol.current.all; stepsol.current.nonmotion = stepsol.current.all;
for i = 1:length(solrb.e.I0)
    stepsol.current.all(:,i) = sum(solrb.e1                          <=solrb.e.I0(i),2)*100/NumNeurons;
    stepsol.opto.all(:,i) = sum(solrb.o1                          <=solrb.o.I0(i),2)*100/NumNeurons;
end

% Options for plotting
options = struct;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.current.all,options); title('Neuron Activation Through Center Electrode'); xlabel('Stimulation Intensity'); ylabel('Percentage Activated Neurons');
ylim([0 100]);
hold on
options.color_area = [128 128 0]./255;
plot_areaerrorbar(stepsol.opto.all,options);
ylim([0 100]);
hold off

%% Half max analysis
e100 = nanmean(solrb.e1);
o100 = nanmean(solrb.o1);
e50 = e100/2;
o50 = o100/2;

Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*e50 + neuron.io.axon(:,ElectrodeNo).*e50; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*o50; % Summation of current directly from stimulus. AU irridance 
frc = Ie_Neurons*.5;
fro = Il_Neurons.*.5;

percent = 0:.05:2;
NumTrials = 100;
numrepeats = 100;
for i = 1:length(percent)
    
    Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*e50.*percent(i) + neuron.io.axon(:,ElectrodeNo).*e50.*percent(i); % Summation of current directly from stimulus + backpropogated up by axons. AU Current
    Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*o50;
    frc = Ie_Neurons*.5;
    fro = Il_Neurons.*.5;
    for ii = 1:numrepeats
        Lambda_Hat_Spikes = Simple_PoissonGen(40+frc+fro, dt, NumTrials);
        Y = prctile(Lambda_Hat_Spikes,05);
        spikee(i,ii) =  Y > 40;
    end
    
    Ie_Neurons = neuron.io.soma(:,ElectrodeNo).*e50 + neuron.io.axon(:,ElectrodeNo).*e50; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
    Il_Neurons = neuron.oo.soma(:,ElectrodeNo).*o50.*percent(i);
    frc = Ie_Neurons*.5;
    fro = Il_Neurons.*.5;
    for ii = 1:numrepeats
        Lambda_Hat_Spikes = Simple_PoissonGen(40+frc+fro, dt, NumTrials);
        Y = prctile(Lambda_Hat_Spikes,05);
        spikeo(i,ii) =  Y > 40;
    end
end
for i=1:length(percent)
    spikee1(i) = sum(spikee(i,:))>50;
    spikeo1(i) = sum(spikeo(i,:))>50;
end

figure;
plot(percent*100,spikee1);
hold on;
plot(percent*100,spikeo1);
xlabel('Percent of Rheobase applied'); ylabel('Spike Detected'); title('Spike 50/50 Determination');
legend('50% Opto','50% Current');
%% Functions
function plot_areaerrorbar(data, options)
    options.color_line = options.color_area;
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