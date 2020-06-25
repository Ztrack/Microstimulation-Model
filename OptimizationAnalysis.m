clc; clear; close all;
load('InitialConditionsFull.mat');

n = 0.5; % Percent activation to count as 'activated'
x = 10; % The current/opto step must activate this many motion tuned neurons to store the step
load('SEoutputc.mat');
threshold.c = thresholdfun(output,units,n,x);
load('SEoutputo.mat');
threshold.o = thresholdfun(output,units,n,x);
%% Microstimulation + Optogenetics Optomization

neuron.inhibitoryfactor = 0.01;
neuron.threshold.rs = 29.0; % Value determined for 95%, used in testing
neuron.threshold.fs = 53.0; % % Value determined for 95%, used in testing
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% lambdatype:
lambdatype = 2; % 3 = MS + Optogenetics (all excitatory - affects all cells indiscriminately)

%% Matlab PSO options
rng default  % For reproducibility
fun = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype);
nvars = 200;
lb = zeros(1,nvars);
ub = [threshold.c threshold.o];
options = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'display','iter');

%% Matlab PSO Once
[x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);

%% Matlab PSO Iterative

options = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'display','off');
numrepeats = 30;
for i = 1:4
    
    lambdatype = i;
    fun = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype);
    for ii = 1:numrepeats
        [x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);
        CostsCurve(i,ii,:) = fval;
    end
    
end

%% Matlab PSO Iterative Across #'s

for i = 1:5:250
    x = i;
    load('SEoutputc.mat');
    threshold.c = thresholdfun(output,units,n,x);
    load('SEoutputo.mat');
    threshold.o = thresholdfun(output,units,n,x);
    for ii = 1:4
        numrepeats = 1;
        for ii = 1:4
            
            lambdatype = ii;
            fun = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype);
            for iii = 1:numrepeats
                [x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);
                CostsCurve(ii,iii,:) = fval;
            end
            
        end
    end
end
        
%% Plotting
CostsCurve(CostsCurve == inf) = NaN;
CostsCurve(CostsCurve == 0) = NaN;
figure;
for i = 1:4
    y = squeeze(nanmean(CostsCurve(i,:,:),2));
    yerr = squeeze(1.96.*nanstd(CostsCurve(i,:,:))/sqrt(numrepeats));
    errorbar(y,yerr);
    hold on;
end
xlabel('Iteration Number');
ylabel('Non-Motion / Motion Neuron Ratio');
title('Optimization Performance');
legend('MS','Opto','MS + Opto (+all)','MS + Opto (-all)');

z = MotionRatioCombined(x,NumNeurons,neuron,lambdatype)
%% Functions
function z = MotionRatioCombined(x,NumNeurons,neuron,lambdatype) % Cost
% z = solution to be minimzied
% ec = all electrode variables
ec = x(1:100); % Electrical stimulation values
eo = x(101:200); % Electrode optical stimulation values

% Electrical / Optical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
Il_Soma_Neurons = zeros(NumNeurons,1); % Initialize luminous intensity

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*neuron.dirmult(:,i);
    Il_Soma_Neurons = Il_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons.
Il_Neurons = Il_Soma_Neurons; % Summation of luminous intensity directly from stimulus.

% Calculate Lambda Hat
[lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,neuron.inhibitoryfactor,lambdatype);

% neuron activation
NumTrials = 100;
simulation = 1;
dt = 1/1000;
bpct = .05;

% for i = 1:NumNeurons
%      output(1,i) = Simple_PoissonGen2(lambdahat(i), dt, NumTrials,simulation,bpct); % Calculate Lambda
% end
% Neuron_Activated = output > neuron.lambda+1;

% TEST METHOD ONLY:
Neuron_Activated = zeros(NumNeurons,1);
Neuron_Activated(neuron.excitatory) = lambdahat(neuron.excitatory) > neuron.threshold.rs;
Neuron_Activated(neuron.inhibitory) = lambdahat(neuron.inhibitory) > neuron.threshold.fs;

% Ratio Calculation

a = sum(Neuron_Activated(neuron.motion.number)); % Number active motion-tuned = activated - inactivated
b = sum(Neuron_Activated(neuron.nonmotion.number)); % Number of active non-motion = activated - inactivated

if a >= length(neuron.motion.number)*.25 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a; % Ratio of what we don't want to activate over what we want to activate
else
    z = inf;
end

end