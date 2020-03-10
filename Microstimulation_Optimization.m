clc; clear; close all;
load('InitialConditionsFull.mat');

%% Microstimulation Optomization

% Problem Variables


% Problem Definiton

problem.CostFunction = @(ec) MotionRatio_MS(ec,NumNeurons,I0_Axon_Neurons,I0_Soma_Neurons,Directional_Current_Mult,Neuron_Motion,Neuron_NonMotion,Neuron_Type,Neuron_Excitatory,Neuron_Inhibitory,Neuron_Motif);  % Cost Function  % Cost Function
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  100;   % Upper Bound of Decision Variables

% Parameters of PSO

params.MaxIt = 30;        % Maximum Number of Iterations
params.nPop = 100;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = 1;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

% Calling PSO

out = PSOFunction(problem, params);

BestSol_MS = out.BestSol;
BestCosts_MS = out.BestCosts;

% Results

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts_MS, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


%% Optogenetics Optimization

% Problem Variables
[~,Neuron_Activated] = MotionRatio_MS(BestSol_MS.Position,NumNeurons,I0_Axon_Neurons,I0_Soma_Neurons,Directional_Current_Mult,Neuron_Motion,Neuron_NonMotion,Neuron_Type,Neuron_Excitatory,Neuron_Inhibitory,Neuron_Motif);

Neuron_No_Activated = find(Neuron_Activated == 1); % All neuron # activated from microstimulation
Neuron_No_Activated_Motion = intersect(Neuron_No_Activated,Neuron_Motion(1,:)); % All neuron # of motion tuned neurons activated
Neuron_No_Activated_NonMotion = intersect(Neuron_No_Activated,Neuron_NonMotion); % All neuron % of non-motion tuned neurons activated
Activated_Motion = length(Neuron_No_Activated_Motion); % Stores # of activated motion-tuned
Activated_NonMotion = length(Neuron_No_Activated_NonMotion); % Stores # of activated nonmotion-tuned
I0_Motion = I0_Soma_Neurons(Neuron_No_Activated_Motion,:); % Stores d information for Motion tuned. sum(I0*e) < 1
I0_NonMotion = I0_Soma_Neurons(Neuron_No_Activated_NonMotion,:); % Stores d information for non-motion tuned. sum(I0*e) > 1
Opto_Thresh = 1; % Threshold that needs to be reached to inactivate cell
% eo = electrode opto stimulus
%%
% Problem Definiton

problem.CostFunction = @(eo) MotionRatio_Opto(eo,Opto_Thresh,I0_Motion,I0_NonMotion,Activated_Motion,Activated_NonMotion);  % Cost Function  % Cost Function
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  100;   % Upper Bound of Decision Variables

% Parameters of PSO

params.MaxIt = 30;        % Maximum Number of Iterations
params.nPop = 100;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

% Calling PSO

out = PSOFunction(problem, params);

BestSol_Opto = out.BestSol;
BestCosts_Opto = out.BestCosts;

% Results

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts_Opto, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


%% Functions
function z = MotionRatio_Opto(eo,Opto_Thresh,I0_Motion,I0_NonMotion,Activated_Motion,Activated_NonMotion)
% z = solution to be minimzied
% eo = all electrode variables

Opto_Motion = zeros(Activated_Motion,1); % Initialize opto stimulus
Opto_NonMotion = zeros(Activated_NonMotion,1); % Initialize opto stimulus

for i = 1:length(eo) % Summation of stimulus for every neuron by every electrode distance & its corresponding stimulus
    Opto_Motion = Opto_Motion + I0_Motion(:,i).*eo(i);
    Opto_NonMotion = Opto_NonMotion + I0_NonMotion(:,i).*eo(i);
end


a = Activated_Motion - sum(Opto_Motion > Opto_Thresh); % Number active motion-tuned = activated - inactivated
b = Activated_NonMotion - sum(Opto_NonMotion > Opto_Thresh); % Number of active non-motion = activated - inactivated

if a >= 40 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a;
else
    z = inf;
end

end

function [z,Neuron_Activated] = MotionRatio_MS(ec,NumNeurons,I0_Axon_Neurons,I0_Soma_Neurons,Directional_Current_Mult,Neuron_Motion,Neuron_NonMotion,Neuron_Type,Neuron_Excitatory,Neuron_Inhibitory,Neuron_Motif)
% z = solution to be minimzied
% ec = all electrode variables
ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
lambda_needed_RS = 29.0; % Value determined experimentally
lambda_needed_FS= 53.0; % % Value determined experimentally

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + I0_Axon_Neurons(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + I0_Soma_Neurons(:,ElectrodeNo(i)).*ec(i).*Directional_Current_Mult(:,ElectrodeNo);
end

Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

% Multiplier Matrices
FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory
Inhibitory_Factor = 0.01;

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection

for i = 1:NumNeurons
    if Neuron_Type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
    end
end

Inhibitory_Effect = zeros(length(Neuron_Inhibitory),1);
for i = 1:length(Neuron_Inhibitory)
    Inhibitory_Effect(i) = Lambda_Hat(Neuron_Inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if Neuron_Type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(Neuron_Motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
    else
        % Inhib Neurons do not inhib themselves
    end
end

Neuron_Activated = zeros(NumNeurons,1);
Neuron_Activated(Neuron_Excitatory) = Lambda_Hat(Neuron_Excitatory) > lambda_needed_RS;
Neuron_Activated(Neuron_Inhibitory) = Lambda_Hat(Neuron_Inhibitory) > lambda_needed_FS;

a = sum(Neuron_Activated(Neuron_Motion(1,:)));
b = sum(Neuron_Activated(Neuron_NonMotion));

if a >= 45 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a; % Minimize # non-motion tuned
else
    z = inf; % If minimum # of motion is not achieved, this is not a good solution
end

end