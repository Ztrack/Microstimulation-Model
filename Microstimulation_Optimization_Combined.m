clc; clear; close all;
load('InitialConditionsFull.mat');

%% Microstimulation + Optogenetics Optomization

neuron.inhibitoryfactor = 0.01;
%neuron.threshold.rs = 29.0; % Value determined for 95%
%neuron.threshold.fs = 53.0; % % Value determined for 95%
neuron.threshold.rs = 21.0; % Value determined experimentally
neuron.threshold.fs = 41.0; % % Value determined experimentally
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% lambdatype:
% 1 = MS + Optogenetics (all silencing),
% 2 = MS + Optogenetics (Inhibitory neuron excitation, only)
% 3 = MS + Optogenetics (Inhibitory neuron excitation, excitatory silencing)
lambdatype = 1;


% Problem Definiton

problem.CostFunction = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype); % Cost
problem.nVar = 200;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  1000;   % Upper Bound of Decision Variables

% Parameters of PSO
params.MaxIt = 100;        % Maximum Number of Iterations
params.AdaptiveItMax = 100; % Max number of adaptive iterations
params.AdaptiveItThreshold = .01; % if the absolute value of it - (it-1) is less than this, we say there is little enough change to end the search
params.nPop = 2000;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = .99;        % Damping Ratio of Inertia Coefficient
params.c1 = 5;              % Personal Acceleration Coefficient
params.c2 = 5;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO - Single

out = PSOFunction(problem, params);

BestSol = out.BestSol;
BestCosts = out.BestCosts;
BestIt = out.NumIt;

figure;
plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Non-Motion / Motion Neuron Ratio');
title('Optimization Performance');
grid on;

%% Calling PSO - Iterative
% it = 1000:1000:1000000;
% for i = 1:length(it)
%     params.nPop = it(i);
%     out = PSOFunction(problem, params);
%     Iterative.Costs(i) = out.BestSol.Cost;
%     BestIt = out.NumIt;
%     
%     if Iterative.Costs(i) == min(Iterative.Costs)
%         Iterative.BestCost = Iterative.Costs(i);
%         Iterative.BestSolPos = out.BestSol.Position;
%     end
% end

% load('Iterative1.mat'); load('Iterative2.mat'); load('Iterative3.mat');
% numit = 20;
% Costs(1,1:numit) = Iterative1.Costs(1:numit);
% Costs(2,1:numit) = Iterative2.Costs(1:numit);
% Costs(3,1:numit) = Iterative3.Costs(1:numit);
% y = mean(Costs);
% ystd = std(Costs);
% err = (ystd./sqrt(size(Costs,1))).*1.96;
% figure; errorbar((1:20)*1000,y,err); title('Number Solutions Affects Optimization Performance'); xlabel('Number Solutions'); ylabel('Optimization Performance');
% 
% figure; plot((1:30)*1000,Iterative3.Costs(1:30));


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
Ir_Soma_Neurons = zeros(NumNeurons,1); % Initialize irridance

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*neuron.dirmult(:,i);
    Ir_Soma_Neurons = Ir_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons. AU Current
Ir_Neurons = Ir_Soma_Neurons; % Summation of current directly from stimulus. AU irridance

% Calculate Lambda Hat
[lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Ir_Neurons,neuron.inhibitoryfactor,lambdatype);

% neuron activation

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