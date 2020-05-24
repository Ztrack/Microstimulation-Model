clc; clear; close all;
load('InitialConditionsFull.mat');

%% Microstimulation Optomization

% Problem Variables


% Problem Definiton

problem.CostFunction = @(ec) MotionRatio_MS(ec,NumNeurons,neuron.io.axon,neuron.io.soma,Directional_Current_Mult,neuron);  % Cost Function  % Cost Function
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  100;   % Upper Bound of Decision Variables

% Parameters of PSO

params.MaxIt = 30;        % Maximum Number of Iterations
params.nPop = 1000;           % Population Size (Swarm Size)
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

figure; set(gcf,'Position',[100 100 800 700]);
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts_MS, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


%% Optogenetics Optimization

% Problem Variables
[~,Neuron_Activated] = MotionRatio_MS(BestSol_MS.Position,NumNeurons,neuron.io.axon,neuron.io.soma,Directional_Current_Mult,neuron);

Neuron_No_Activated = find(Neuron_Activated == 1); % All neuron # activated from microstimulation
Neuron_No_Activated_Motion = intersect(Neuron_No_Activated,neuron.motion.number); % All neuron # of motion tuned neurons activated
Neuron_No_Activated_NonMotion = intersect(Neuron_No_Activated,neuron.nonMotion.number); % All neuron % of non-motion tuned neurons activated
Activated_Motion = length(Neuron_No_Activated_Motion); % Stores # of activated motion-tuned
Activated_NonMotion = length(Neuron_No_Activated_NonMotion); % Stores # of activated nonmotion-tuned
I0_Motion = neuron.io.soma(Neuron_No_Activated_Motion,:); % Stores d information for Motion tuned. sum(I0*e) < 1
I0_NonMotion = neuron.io.soma(Neuron_No_Activated_NonMotion,:); % Stores d information for non-motion tuned. sum(I0*e) > 1
Opto_Thresh = 1; % Threshold that needs to be reached to inactivate cell
% eo = electrode opto stimulus

%%
% Problem Definiton

problem.CostFunction = @(eo) MotionRatio_Opto(eo,Opto_Thresh,I0_Motion,I0_NonMotion,Activated_Motion,Activated_NonMotion,neuron);  % Cost Function  % Cost Function
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  100;   % Upper Bound of Decision Variables

% Parameters of PSO

params.MaxIt = 30;        % Maximum Number of Iterations
params.nPop = 100000;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = 1;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

% Calling PSO

out = PSOFunction(problem, params);

BestSol_Opto = out.BestSol;
BestCosts_Opto = out.BestCosts;

% Results

figure; set(gcf,'Position',[100 100 800 700]);
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts_Opto, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

% save('BestSol.mat','BestSol_MS','BestCosts_MS','BestSol_Opto','BestCosts_Opto');
% load('BestSol.mat');
%% Stimulus Plots

Stim_Current = zeros(sx, sy);
Stim_Opto = zeros(sx, sy);
for i = 1:length(electrode.x)
    Stim_Loc =  zeros(sx, sy); % Location of stimulus
    Stim_Loc(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1; % Builds electrode stimulus location matrix
    Ed = (bwdist(Stim_Loc)); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Ed(Ed == 0) = 1;
    R2 = 1./(Ed);
    % R2(R2 < 0.01) = 0; % Filtering
    Stim_Current = Stim_Current + R2.*BestSol_MS.Position(i); % Distance of every square to the nearest stim location
    Stim_Opto = Stim_Opto + R2.*BestSol_Opto.Position(i);
end
grad=colorGradient([1/4 1/2 1],[0 0 1],128); map = [1 1 1; grad];
figure; set(gcf,'Position',[100 100 800 700]); imagesc(Stim_Current); colorbar; title('Current Intensity'); colormap(map);
figure; set(gcf,'Position',[100 100 800 700]); imagesc(Stim_Opto); colorbar; title('Light Intensity'); colormap(map);

%% Neuron Pop Plots
[Lambda,Lambda_Hat] = Lambda_Hat_Gen(BestSol_MS,BestSol_Opto,NumNeurons,Directional_Current_Mult,neuron);
Neuron_Pop = zeros(sx, sy); Neuron_Pop_LH = Neuron_Pop; Neuron_Pop_L = Neuron_Pop;
for i = 1:NumNeurons
    Neuron_Pop(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii,neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = 1;
    Neuron_Pop_LH(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii,neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = Lambda_Hat(i);
    Neuron_Pop_L(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii,neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = Lambda(i);
end
Neuron_Pop_D = Neuron_Pop_LH - Neuron_Pop_L;
Neuron_Pop_D(Neuron_Pop_D < 0) = 0; % Neuron pop delta firing rate

figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_D); colorbar; title('Neural Population Firing Rate Change'); colormap(map);

%%
Stim_Current1 = Stim_Current;
Stim_Current1(Stim_Current1 > 1) = 1;
Stim_Opto1 = Stim_Opto;
Stim_Opto1(Stim_Opto1 > 1) = 1;

[y1,x1] = ind2sub([sx sy],population.inhibitory.indices); % Inhibitory
[y2,x2] = ind2sub([sx sy],population.excitatory.indices); % Non-Motion Excitatory
[y3,x3] = ind2sub([sx sy],population.motion.indices); % Motion
[y4,x4] = ind2sub([sx sy],population.axon.indices); % Axons

grad=colorGradient([1 1 1],[0 0 1],128); map = [1 1 1; grad];

figure; set(gcf,'Position',[100 100 800 700]); % Stim Current with neuron figure
imagesc(Stim_Current1); colormap(map);
hold on
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Current Stimulus Over Neural Population'); %legend('Inhibitory','Excitatory');

figure; set(gcf,'Position',[100 100 800 700]); % Stim Opto with neuron figure
imagesc(Stim_Opto1); colormap(map);
hold on
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Light Stimulus Over Neural Population'); %legend('Inhibitory','Excitatory');

%% Functions

function [z,Neuron_Activated] = MotionRatio_MS(ec,NumNeurons,Directional_Current_Mult,neuron)
% z = solution to be minimzied
% ec = all electrode variables
ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
lambda_needed_RS = 29.0; % Value determined experimentally
lambda_needed_FS= 53.0; % % Value determined experimentally

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*Directional_Current_Mult(:,ElectrodeNo);
end

Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

% Multiplier Matrices
FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory
Inhibitory_Factor = 0.01;

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection

for i = 1:NumNeurons
    if neuron.type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
    end
end

Inhibitory_Effect = zeros(length(neuron.inhibitory),1);
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(i) = Lambda_Hat(neuron.inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if neuron.type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(neuron.motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
    else
        % Inhib Neurons do not inhib themselves
    end
end

Neuron_Activated = zeros(NumNeurons,1);
Neuron_Activated(neuron.excitatory) = Lambda_Hat(neuron.excitatory) > lambda_needed_RS;
Neuron_Activated(neuron.inhibitory) = Lambda_Hat(neuron.inhibitory) > lambda_needed_FS;

a = sum(Neuron_Activated(neuron.motion.number));
b = sum(Neuron_Activated(neuron.nonMotion.number));

if a >= length(neuron.motion.number)*.25 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a; % Minimize # non-motion tuned
else
    z = inf; % If minimum # of motion is not achieved, this is not a good solution
end

end

function z = MotionRatio_Opto(eo,Opto_Thresh,I0_Motion,I0_NonMotion,Activated_Motion,Activated_NonMotion,neuron)
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

if a >= length(neuron.motion.number)*.20 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a;
else
    z = inf;
end

end

function [Lambda,Lambda_Hat] = Lambda_Hat_Gen(BestSol_MS,BestSol_Opto,NumNeurons,Directional_Current_Mult,neuron)

% ec = all electrode current variables
% eo = all electrode light variables
ec = BestSol_MS.Position;
eo = BestSol_Opto.Position;

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
lambda_needed_RS = 29.0; % Value determined experimentally
lambda_needed_FS= 53.0; % % Value determined experimentally

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*Directional_Current_Mult(:,ElectrodeNo);
end

Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

% Multiplier Matrices
FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory
Inhibitory_Factor = 0.01;

Lambda = zeros(NumNeurons,1);
Lambda(neuron.excitatory) = RS_Lambda;
lambda(neuron.inhibitory) = FS_Lambda;

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection

for i = 1:NumNeurons
    if neuron.type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
    end
end

Inhibitory_Effect = zeros(length(neuron.inhibitory),1);
for i = 1:length(neuron.inhibitory)
    Inhibitory_Effect(i) = Lambda_Hat(neuron.inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if neuron.type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(neuron.motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
    else
        % Inhib Neurons do not inhib themselves
    end
end

Opto_Thresh = 1;
Opto_Stim = zeros(NumNeurons,1); % Initialize opto stimulus

for i = 1:length(eo) % Summation of stimulus for every neuron by every electrode distance & its corresponding stimulus
    Opto_Stim = Opto_Stim + neuron.io.soma(:,i).*eo(i);
end

for i = 1:NumNeurons
    if Opto_Stim(i) >= Opto_Thresh
        Lambda_Hat(i) = 0;
    end
end

Lambda_Hat(Lambda_Hat > 500) = 500;
Lambda_Hat(Lambda_Hat < 0) = 0;
end

function [grad,im]=colorGradient(c1,c2,depth)
% COLORGRADIENT allows you to generate a gradient between 2 given colors,
% that can be used as colormap in your figures.
%
% USAGE:
%
% [grad,im]=getGradient(c1,c2,depth)
%
% INPUT:
% - c1: color vector given as Intensity or RGB color. Initial value.
% - c2: same as c1. This is the final value of the gradient.
% - depth: number of colors or elements of the gradient.
%
% OUTPUT:
% - grad: a matrix of depth*3 elements containing colormap (or gradient).
% - im: a depth*20*3 RGB image that can be used to display the result.
%
% EXAMPLES:
% grad=colorGradient([1 0 0],[0.5 0.8 1],128);
% surf(peaks)
% colormap(grad);
%
% --------------------
% [grad,im]=colorGradient([1 0 0],[0.5 0.8 1],128);
% image(im); %display an image with the color gradient.
% Copyright 2011. Jose Maria Garcia-Valdecasas Bernal
% v:1.0 22 May 2011. Initial release.
%Check input arguments.
%input arguments must be 2 or 3.
error(nargchk(2, 3, nargin));
%If c1 or c2 is not a valid RGB vector return an error.
if numel(c1)~=3
    error('color c1 is not a valir RGB vector');
end
if numel(c2)~=3
    error('color c2 is not a valir RGB vector');
end
if max(c1)>1&&max(c1)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c1 is not given as intensity values. Trying to convert');
    c1=c1./255;
elseif max(c1)>255||min(c1)<0
    error('C1 RGB values are not valid.')
end
if max(c2)>1&&max(c2)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c2 is not given as intensity values. Trying to convert');
    c2=c2./255;
elseif max(c2)>255||min(c2)<0
    error('C2 RGB values are not valid.')
end
%default depth is 64 colors. Just in case we did not define that argument.
if nargin < 3
    depth=64;
end
%determine increment step for each color channel.
dr=(c2(1)-c1(1))/(depth-1);
dg=(c2(2)-c1(2))/(depth-1);
db=(c2(3)-c1(3))/(depth-1);
%initialize gradient matrix.
grad=zeros(depth,3);
%initialize matrix for each color. Needed for the image. Size 20*depth.
r=zeros(20,depth);
g=zeros(20,depth);
b=zeros(20,depth);
%for each color step, increase/reduce the value of Intensity data.
for j=1:depth
    grad(j,1)=c1(1)+dr*(j-1);
    grad(j,2)=c1(2)+dg*(j-1);
    grad(j,3)=c1(3)+db*(j-1);
    r(:,j)=grad(j,1);
    g(:,j)=grad(j,2);
    b(:,j)=grad(j,3);
end
%merge R G B matrix and obtain our image.
im=cat(3,r,g,b);
end