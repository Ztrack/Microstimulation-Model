clearvars; clc;
load('InitialConditionsFullB.mat');
% This initial condition has 5 Excitatory and 5 Inhibitory

%%

% Features
calctype = 2;
Motion_Axon = 0; % if set =1, enabled, connections between relevant motion neurons is established
Oscillatory_Behavior = 1; % is set =1, .33 of all inhibitory neurons will use Oscillatory poisson function
Directional_Current_Modifier = 1; % if set =1 & enabled, multiplier is applied to the soma depending on axon-hillock location
lambdatype = 1;

% Apply Features
if Directional_Current_Modifier == 0
    neuron.dirmult(:,:) = 1;
end
if Motion_Axon == 0
    I0_Motion_Neurons = zeros(NumNeurons,length(electrode.x));
end

% Parameters

bpct = 50; % Bottom Percentile for Rheobase calculation. 50 for 50%, 05 for 95% CI.
inhibitoryfactor = [0.01 0.005 0.001 0.0005 0.0001];
NumTrials = 100;

%% Loop Start
h = 100; % Step size
I0 = 0:h:100000;  % Current Steps
numrepeats = 12; % Number of overall repeats
ElectrodeDist = sqrt((sx/2-electrode.x).^2 + (sy/2-electrode.y).^2);
ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode
 
for kkk = 1:2
    lambdatype = kkk;
for kk = 1:length(inhibitoryfactor)
    load('InitialConditionsFullB.mat');
    
    NumInhibitoryMotif = 2; % Number of inhibitory neurons
    NumNeuronsMotif = 5+NumInhibitoryMotif; % Reinitializes # Neurons per motif
    NumNeurons = NumMotifs*NumNeuronsMotif;
    a = zeros(1,NumInhibitoryMotif); b= ones(1,5-NumInhibitoryMotif); a = [a b]; a = [a 0 0 0 0 0];
    Extra_Inhib_Neurons = repmat(a,1,NumMotifs); Extra_Inhib_Neurons = find(Extra_Inhib_Neurons == 1); % Identifies extra neurons
    Neuron_RB = NaN(numrepeats,NumNeurons); % Rhoebase for every neuron, stored as I0 which causes neuron.lambda+1 spike
    neuron.io.axon(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
    neuron.io.soma(Extra_Inhib_Neurons,:) = []; % Erases data from extra neurons
    neuron.motif(Extra_Inhib_Neurons) = [];
    neuron.type(Extra_Inhib_Neurons) = [];
    neuron.oo.soma(Extra_Inhib_Neurons,:) = [];
    neuron.dirmult(Extra_Inhib_Neurons,:) = [];
    neuron.lambda = zeros(1,NumNeurons);
    neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
    neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
    neuron.inhibitory = find(neuron.type == 1); % New Inhibitory List
    neuron.excitatory = find(neuron.type == 2);
    parfor jj = 1:numrepeats
        
        Neuron_RB1 = NaN(1,NumNeurons);
        
        for ii = 1:length(I0)
            Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
            Ie_Soma_Neurons = zeros(NumNeurons,1);
            Ir_Soma_Neurons = zeros(NumNeurons,1); % Initialize irridance
            
            for i = 1:length(ElectrodeNo) % Summation of current for every neuron component by every electrode & its corresponding current
                Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i));
                Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*neuron.dirmult(:,i);
                Ir_Soma_Neurons = Ir_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i));
            end

            Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
            Ir_Neurons = Ir_Soma_Neurons; % Summation of current directly from stimulus. AU irridance
            
            % Calculate neuron.lambda Change for poisson process
            [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Ir_Neurons,inhibitoryfactor(kk),lambdatype);
            
            % Finding RB for each neuron
            
            for i = 1:NumNeurons
                if isnan(Neuron_RB1(i)) % If RB does not exist, continue, otherwise this neuron is skipped
                    if neuron.oscillatorytype(i) == 1 & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
                        Lambda_Hat_Spikes = Oscillatory_PoissonGen(lambdahat(i), dt, NumTrials);
                    else % Otherwise use the simple function:
                        Lambda_Hat_Spikes = Simple_PoissonGen(lambdahat(i), dt, NumTrials);
                    end
                    
                    Y = prctile(Lambda_Hat_Spikes,bpct); % Calculates bottom xth percentile
                    if Y > mean(neuron.lambda(i)+1)
                        Neuron_RB1(i) = I0(ii);
                    end
                end
            end
        end
        Neuron_RB(jj,:) = Neuron_RB1;
    end
    
    %% Export
    if lambdatype == 1
        if kk == 1
            solrbb.e1 = Neuron_RB;
        elseif kk ==2
            solrbb.e2 = Neuron_RB;
        elseif kk ==3
            solrbb.e3 = Neuron_RB;
        elseif kk ==4
            solrbb.e4 = Neuron_RB;
        elseif kk ==5
            solrbb.e5 = Neuron_RB;
        elseif kk ==6
            solrbb.e5 = Neuron_RB;
        end
    elseif lambdatype == 2
        if kk == 1
            solrbb.o1 = Neuron_RB;
        elseif kk ==2
            solrbb.o2 = Neuron_RB;
        elseif kk ==3
            solrbb.o3 = Neuron_RB;
        elseif kk ==4
            solrbb.o4 = Neuron_RB;
        elseif kk ==5
            solrbb.o5 = Neuron_RB;
        end
    end
end
end
if calctype == 1
    save('solrb1.mat','solrb','I0','h','numrepeats','ElectrodeNo');
else
    save('solrbb.mat','solrb','I0','h','numrepeats','ElectrodeNo');
end