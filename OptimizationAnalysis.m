clc; clear; close all;
load('InitialConditionsFull.mat');
load('threshold.mat');
%% Microstimulation + Optogenetics Optomization

neuron.inhibitoryfactor = 0.01;
neuron.threshold.rs = 29.0; % Value determined for 95%
neuron.threshold.fs = 53.0; % % Value determined for 95%
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% lambdatype:
lambdatype = 7; % 3 = MS + Optogenetics (all excitatory - affects all cells indiscriminately)


%% Problem Definiton

problem.CostFunction = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype); % Cost
problem.nVar = 200;       % Number of Unknown (Decision) Variables
problem.VarMin =  zeros(200,1);
problem.VarMax =  [threshold.c threshold.o]';   % Upper Bound of Decision Variables

% Parameters of PSO
params.MaxIt = 50;        % Maximum Number of Iterations
params.AdaptiveItMax = 50; % Max number of adaptive iterations
params.AdaptiveItThreshold = .01; % if the absolute value of it - (it-1) is less than this, we say there is little enough change to end the search
params.nPop = 100;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = .99;        % Damping Ratio of Inertia Coefficient
params.c1 = 5;              % Personal Acceleration Coefficient
params.c2 = 5;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin


%% Calling PSO
out = PSOFunction2(problem, params);


%% Calling PSO

numrepeats = 100;
BestCost = inf;
for i = 1:7
    
    lambdatype = i;
    for ii = 1:numrepeats
        out = PSOFunction2(problem, params);
        CostsCurve(i,ii,:) = out.BestCosts;
    end
    
end

% figure;
% plot(CostsCurve, 'LineWidth', 2);
% semilogy(CostsCurve, 'LineWidth', 2);
% xlabel('Iteration Number');
% ylabel('Non-Motion / Motion Neuron Ratio');
% title('Optimization Performance');
% grid on;

save optimize.mat;
%% Functions
function z = MotionRatioCombined(x,NumNeurons,neuron,lambdatype) % Cost
% z = solution to be minimzied
% ec = all electrode variables
ec = x(1:100); % Electrical stimulation values
eo = x(101:200); % Electrode optical stimulation values

% Electrical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
Il_Soma_Neurons = zeros(NumNeurons,1); % Initialize irridance

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*neuron.dirmult(:,i);
    Il_Soma_Neurons = Il_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
Il_Neurons = Il_Soma_Neurons; % Summation of current directly from stimulus. AU irridance

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