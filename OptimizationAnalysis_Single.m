clc; clear; close all;
load('InitialConditionsFull.mat');

%% Microstimulation Optomization

% Problem Variables
neuron.inhibitoryfactor = 0.01;
neuron.threshold.rs = 29.0; % Value determined for 95%, used in testing
neuron.threshold.fs = 53.0; % % Value determined for 95%, used in testing
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% lambdatype:
lambdatype = 1; % 3 = MS + Optogenetics (all excitatory - affects all cells indiscriminately)

% Problem Definiton

problem.CostFunction = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype); % Cost
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  1000;   % Upper Bound of Decision Variables

% Parameters of PSO

params.MaxIt = 30;        % Maximum Number of Iterations
params.AdaptiveItMax = 20;
params.AdaptiveItThreshold = .01;
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

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts_MS, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


%% Calculate lambdahat for above

ec = BestSol_MS.Position; % Electrical stimulation values

% Electrical / Optical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
Il_Soma_Neurons = zeros(NumNeurons,1); % Initialize luminous intensity

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*neuron.dirmult(:,i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons.
Il_Neurons = Il_Soma_Neurons; % Summation of luminous intensity directly from stimulus.

% Calculate Lambda Hat
[lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,neuron.inhibitoryfactor,lambdatype);
neuron.mslambdahat = lambdahat;

%% Optogenetics Optimization

% Problem Definiton
lambdatype = 4;

problem.CostFunction = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype); % Cost
problem.nVar = 100;       % Number of Unknown (Decision) Variables
problem.VarMin =  0;  % Lower Bound of Decision Variables
problem.VarMax =  1;   % Upper Bound of Decision Variables

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

function z = MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype) % Cost
% z = solution to be minimzied
% ec = all electrode variables
ec = electrodestim(1:100); % Electrical stimulation values
eo = zeros(100,1);


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

if lambdatype == 3
    [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,neuron.inhibitoryfactor,2);
    lambdahat = neuron.mslambdahat + (lambdahat - neuron.lambda);
elseif lambdatype == 4
    [lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,neuron.inhibitoryfactor,2);
    lambdahat = neuron.mslambdahat - (lambdahat - neuron.lambda);
end

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