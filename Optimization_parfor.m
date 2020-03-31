clc; clear; close all;
load('InitialConditionsFull.mat');

%% Microstimulation + Optogenetics Optomization

% Problem Definiton

problem.CostFunction = @(x) MotionRatioCombined(x,NumNeurons,neuron,Directional_Current_Mult); % Cost
problem.nVar = 200;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  50;   % Upper Bound of Decision Variables

% Parameters of PSO
params.AdaptiveItMax = 5; % Max number of adaptive iterations
params.AdaptiveItThreshold = .01; % if the absolute value of it - (it-1) is less than this, we say there is little enough change to end the search
params.MaxIt = 50;        % Maximum Number of Iterations
params.nPop = 1000;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = 1;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

% Calling PSO Iterative

it = 1000:1000:20000;
parfor ii = 1:50
    Iterative_Costs1 = [];
    Iterative_BestCost1 = [];
    Iterative_BestSolPos1 = [];
for i = 1:size(it,2)
    nPop = it(i);
    out = PSOFunctionparfor(problem, params,nPop);
    Iterative_Costs1 = out.BestSol.Cost;
    
    if Iterative_Costs1(i) == min(Iterative_Costs1)
        Iterative_BestCost1 = Iterative_Costs1(i);
        Iterative_BestSolPos1 = out.BestSol.Position;
    end
end
Iterative_Costs(ii,:) = Iterative_Costs1;
Iterative_BestCost(ii,:) = Iterative_BestCost1;
Iterative_BestSolPos(ii,:) = Iterative_BestSolPos1;
end

% Iterative1 = Iterative
% save('Iterative','Iterative');
% Plotting:
% figure;
% % plot(BestCosts, 'LineWidth', 2);
% semilogy(BestCosts, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;

%%
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

%%

%% Functions
function z = MotionRatioCombined(x,NumNeurons,neuron,Directional_Current_Mult) % Cost
% z = solution to be minimzied
% ec = all electrode variables
ec = x(1:100);
eo = x(101:200);

% Electrical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
Ir_Soma_Neurons = zeros(NumNeurons,1); % Initialize irridance

lambda_needed_RS = 29.0; % Value determined experimentally
lambda_needed_FS= 53.0; % % Value determined experimentally

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*Directional_Current_Mult(:,ElectrodeNo);
    Ir_Soma_Neurons = Ir_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons
Ir_Neurons = Ir_Soma_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

% Calculate Lambda Hat

FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory
Inhibitory_Factor = 0.01;

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection

for i = 1:NumNeurons
    if neuron.type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Neurons(i) - (FS_Lambda * Ir_Neurons(i))); % Lambda_Hat firing rate based off microstimulation
    end
end

Inhibitory_Effect = zeros(length(neuron.inhibitory),1);
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(i) = Lambda_Hat(neuron.inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if neuron.type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Neurons(i))) - (FS_Lambda * Ir_Neurons(i)) - Inhibitory_Effect(neuron.motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
    else
        % Inhib Neurons do not inhib themselves
    end
end

% neuron activation

Neuron_Activated = zeros(NumNeurons,1);
Neuron_Activated(neuron.excitatory) = Lambda_Hat(neuron.excitatory) > lambda_needed_RS;
Neuron_Activated(neuron.inhibitory) = Lambda_Hat(neuron.inhibitory) > lambda_needed_FS;

% Ratio Calculation

a = sum(Neuron_Activated(neuron.motion.number)); % Number active motion-tuned = activated - inactivated
b = sum(Neuron_Activated(neuron.nonmotion.number)); % Number of active non-motion = activated - inactivated

if a >= 45 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a;
else
    z = inf;
end

end