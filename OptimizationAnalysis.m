clc; clear; close all;
load('InitialConditionsFull.mat');
load('threshold.mat');
%% Microstimulation + Optogenetics Optomization

neuron.inhibitoryfactor = 0.01;
%neuron.threshold.rs = 29.0; % Value determined for 95%
%neuron.threshold.fs = 53.0; % % Value determined for 95%
neuron.threshold.rs = 21.0; % Value determined experimentally
neuron.threshold.fs = 41.0; % % Value determined experimentally
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

h = 50; % number of steps for current/LI
unitsmax(1) = 30000; % Point at which 100% of neurons are activated
units50(1) = 10000; % Point at which 50% of neurons are activated
unitsmax(2) = 100000;
units50(2) = 30000; 

u1 = linspace(0,units50(1),h*.8);  % Current OR liminous intensity Steps
u2 = linspace(0,units50(2),h*.8);
u1 = [u1 linspace(units50(1)+unitsmax(1)*.2*h,unitsmax(1),h*.2)];
u2 = [u2 linspace(units50(2)+unitsmax(2)*.2*h,unitsmax(2),h*.2)];

% lambdatype:
lambdatype = 1; % 3 = MS + Optogenetics (all excitatory - affects all cells indiscriminately)


%% Problem Definiton

problem.CostFunction = @(x) MotionRatioCombined(x,NumNeurons,neuron,lambdatype); % Cost
problem.nVar = 200;       % Number of Unknown (Decision) Variables
problem.VarMin =  [threshold.c-1.*1250 threshold.o-1.*1500];
problem.VarMin(problem.VarMin < 0) = 0; % Fix non-zeros
problem.VarMax =  [threshold.c .*1250 threshold.o .*100];   % Upper Bound of Decision Variables

%problem.VarMin =  [0:99 15000:15099];  % TEST VARIABLES - 0
%problem.VarMax =  [1000:1099 15000:15099]; % TEST VARIABLES

% Parameters of PSO
params.MaxIt = 50;        % Maximum Number of Iterations
params.AdaptiveItMax = 50; % Max number of adaptive iterations
params.AdaptiveItThreshold = .01; % if the absolute value of it - (it-1) is less than this, we say there is little enough change to end the search
params.nPop = 2000;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = .99;        % Damping Ratio of Inertia Coefficient
params.c1 = 5;              % Personal Acceleration Coefficient
params.c2 = 5;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin


%% Calling PSO

numrepeats = 100;
BestCost = inf;

for ii = 1:numrepeats
    out = PSOFunction2(problem, params);
    
    if out.BestSol.Cost < BestCost
        BestCost = out.BestSol.Cost;
        BestPos = out.BestSol.Position;
        CostsCurve = out.BestCosts;
    end
end

figure;
plot(CostsCurve, 'LineWidth', 2);
semilogy(CostsCurve, 'LineWidth', 2);
xlabel('Iteration Number');
ylabel('Non-Motion / Motion Neuron Ratio');
title('Optimization Performance');
grid on;

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