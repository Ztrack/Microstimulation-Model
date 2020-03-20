clearvars; clc; close all;
load('InitialConditionsFullB.mat')
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)

Stim_Center = 0; % If 1, uses a stimulus at center electrode instead of random
Num_Stim = 5;
I0 = 100;

%% Stimulation Locations

if Stim_Center == 1
    ElectrodeDist = sqrt((sx/2-electrode.x).^2 + (sy/2-electrode.y).^2);
    ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode
else
    ElectrodeNo = randperm(100,Num_Stim);
end

Stim_Loc =  zeros(sx, sy);
for i = 1:length(ElectrodeNo)
Stim_Loc(electrode.y(ElectrodeNo(i))-electrode.radii:electrode.y(ElectrodeNo(i))+electrode.radii,electrode.x(ElectrodeNo(i))-electrode.radii:electrode.x(ElectrodeNo(i))+electrode.radii)  = 1;
end
Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
%% Extracellular Current Matrix


Ie_Axon_Neurons = sum(neuron.io.axon(:,ElectrodeNo).*I0,2);
Ie_Soma_Neurons = sum(neuron.io.soma(:,ElectrodeNo).*I0,2);
Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons

%% Multiplier Matrices 
% Constructs effect of microstimulation on neuron firing
dt = 1/1000; % step/bin time length 1ms
t_total= 0:dt:T-dt ;
NumSteps = 1/dt; % Number of steps/bins in time t_total
FS_Lambda = 40; % Fast Spiking Neurons, Inhibitory
RS_Lambda = 20; % Regular Spiking Neurons, Excitory

Lambda_Hat = zeros(NumNeurons,1); % Probability change due to current injection
Lambda = zeros(NumNeurons,1);

for i = 1:NumNeurons
    if neuron.type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
        Lambda(i) = FS_Lambda; % Lambda regular, no microstimulation
    end
end

Inhibitory_Effect = zeros(NumMotifs,1);
Inhibitory_Effect_Lambda = zeros(NumMotifs,1);
for i = 1:NumMotifs
    Inhibitory_Effect(i) = Lambda_Hat(neuron.inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
    Inhibitory_Effect_Lambda = FS_Lambda*Inhibitory_Factor;
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if neuron.type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(neuron.motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
        Lambda(i) = RS_Lambda - Inhibitory_Effect_Lambda; % Lambda regular, no microstimulation
    else
        % Inhib Neurons do not inhib themselves
    end
end

%% Poisson Spiking Model pt 1 - Lambda (No Stimulation)
Distances = zeros(length(ElectrodeNo),length(motif.center.x));
for i = 1:length(ElectrodeNo)
    Distances(i,:) = sqrt((electrode.x(ElectrodeNo(i))-motif.center.x).^2 + (electrode.y(ElectrodeNo(i)) - motif.center.y).^2); % Stores distances of every motif to every active electrode
end
Dist = sort(Distances(1,:));
Stimulated_Motifs = find(Distances(1,:) == Dist(1)); % Stores which motif is closest
Stimulated_Neurons = find(neuron.motif == Stimulated_Motifs(1));
    
NumTrials = 1000;
Oscillatory_Behavior = 1;
for i = 1:NumNeurons
    if intersect(neuron.oscillatory,i) == i & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
        Lambda_Spikes = Oscillatory_PoissonGen(Lambda(i), dt, NumTrials);
    else % Otherwise use the simple function:
        Lambda_Spikes = Simple_PoissonGen(Lambda(i), dt, NumTrials);
    end
    Lambda_Hat_Spikes_Mean(i) =  mean(Lambda_Spikes);
end

Neuron_Pop_Spikes_Lambda = zeros(sx, sy);
for i = 1:NumNeurons
    Neuron_Pop_Spikes_Lambda(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii, neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = Lambda_Spikes(i);
end

%% Poisson Spiking Model pt 2 - Lambda Hat (Stimulated)
NumTrials = 1000;
Oscillatory_Behavior = 1;
for i = 1:NumNeurons
    if intersect(neuron.oscillatory,i) == i & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
        Lambda_Hat_Spikes = Oscillatory_PoissonGen(Lambda_Hat(i), dt, NumTrials);
    else % Otherwise use the simple function:
        Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
    end
    Lambda_Hat_Spikes_Mean(i) =  mean(Lambda_Hat_Spikes);
end


Neuron_Pop_Spikes = zeros(sx, sy);
for i = 1:length(neuron.x)
    Neuron_Pop_Spikes(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii, neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = Lambda_Hat_Spikes_Mean(i);
end

%% Plots for Model

figure; set(gcf,'Position',[100 100 800 700]);
map = [1 1 1; 0 0 0]; colormap(map);
imagesc(1./Stim_Distance); title('Active Electrode Locations'); xlabel('X Position (mm)'); ylabel('Y Position (mm)');

[y1,x1] = ind2sub([sx sy],population.inhibitory.indices); % Inhibitory
[y2,x2] = ind2sub([sx sy],population.excitatory.indices); % Non-Motion Excitatory
[y3,x3] = ind2sub([sx sy],population.motion.indices); % Motion
[y4,x4] = ind2sub([sx sy],population.axon.indices); % Axons

population.empty = zeros(sx,sy);
figure; set(gcf,'Position',[100 100 800 700]); 
map = [1 1 1];
imagesc(population.empty); colormap(map);
hold on
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); %legend('Inhibitory','Excitatory');

%% Colormaps

% Colormap creation
map = [0 0 0; 239 122 37; 81 182 106; 50 127 183; 152 85 159; 225 152 153;]./255; % Colormap. First number is 0, then 1 to 15.
c1 = .85;
c2 = .7;
map = [map ; map(2,1:3)*c1 ; map(3,1:3)*c1 ; map(4,1:3)*c1 ; map(5,1:3)*c1 ; map(6,1:3)*c1];
map = [map ; 1 1 1 ; map(3,1:3)*c2 ; map(4,1:3)*c2 ; map(5,1:3)*c2 ; map(6,1:3)*c2];

population.pads.map = zeros(sx,sy);
for i = 1:NumPads
population.pads.map(pad.start.y(i):pad.end.y(i),pad.start.x(i):pad.end.x(i)) = i;
end

figure; set(gcf,'Position',[100 100 800 700]); imagesc(population.pads.map); colormap(map); title('Pads'); % Color coded map representing different pads
figure; set(gcf,'Position',[100 100 800 700]); imagesc(population.pads.map); colormap(map); title('Pads with electrodes'); % Color coded map representing different pads
hold on;
plot(electrode.x,electrode.y,'.','color','black');

% Neuron Connections

% figure; set(gcf,'Position',[100 100 800 700]); imagesc(population.pads.map); colormap(map); title('Neuron Synaptic Connections'); % Color coded map representing different pads
% hold on
% for i = 1:length(Neuron_Connected)
%         x1 = neuron.x(Neuron_Connected(i,1)); x2 = neuron.x(Neuron_Connected(i,2));
%         y1 = neuron.y(Neuron_Connected(i,1)); y2 = neuron.y(Neuron_Connected(i,2));
%         line([x1,x2], [y1,y2],'Color', 'black');
% end

%%
grad=colorGradient([0 0 1],[1 1 0],128); map = [1 1 1; grad];

% Stimulated motifs / neuron pop maps:
Stimulated_Neurons = find(Lambda_Hat > 50 & Lambda_Hat < 100);
Stimulated_Motifs = neuron.motif(Stimulated_Neurons);

figure; set(gcf,'Position',[100 100 800 700]); 
imagesc(Neuron_Pop_Spikes); 
title('Stimulated Motif'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 max(Lambda_Hat(Stimulated_Neurons))]); colormap(map);
xlim([motif.center.x(Stimulated_Motifs(1))-50 motif.center.x(Stimulated_Motifs(1))+50]); ylim([motif.center.y(Stimulated_Motifs(1))-50 motif.center.y(Stimulated_Motifs(1))+50]);
hold on
plot(x4,y4,'.','color','Black');


%%
figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_Spikes); title('Neuron Response Stimulated'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 100]);
colormap(map);

% Non-Stimulated motifs / neuron pop maps:
figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_Spikes_Lambda(rangey,rangex)); title('Non-Stimulated Motif'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 max(Lambda_Hat(Stimulated_Neurons))]);
hold on
a = find(population.empty(rangey,rangex) > 0);
[y1,x1] = ind2sub([length(rangex) length(rangey)],a);
plot(x1,y1,'.','color','Black');
colormap(map);

figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_Spikes_Lambda); title('Neuron Response Non-Stimulated'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar;
colormap(map);

% motifs difference map:
figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_Spikes(rangey,rangex) - Neuron_Pop_Spikes_Lambda(rangey,rangex)); title('Motif Difference Map'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([0 max(Lambda_Hat(Stimulated_Neurons)-max(Lambda(Stimulated_Neurons)))]);
hold on
a = find(population.empty(rangey,rangex) > 0);
[y1,x1] = ind2sub([length(rangex) length(rangey)],a);
plot(x1,y1,'.','color','Black')
colormap(map);

figure; set(gcf,'Position',[100 100 800 700]); imagesc(Neuron_Pop_Spikes-Neuron_Pop_Spikes_Lambda); title('Neuron Response Difference'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 100]);
colormap(map);

%% Functions

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

function Trial_Spikes = Simple_PoissonGen(Lambda, dt, NumTrials)
NumSteps = 1/dt;
Spike_Probability = Lambda.*dt;
X_random = rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end

function Trial_Spikes = Oscillatory_PoissonGen(freqOscillation, dt, NumTrials)
dur = 1;
modIndex = 1;
NumSteps = 1/dt;
t = 1:1:NumSteps;
Spike_Probability = (freqOscillation * (1/dur) * (1 + modIndex .* cos(2*pi * freqOscillation * t/NumSteps))) * dt;
X_random = 1 + (0-1).*rand(NumTrials,NumSteps);
Spikes = (X_random < Spike_Probability);
Trial_Spikes = sum(Spikes, 2);
end