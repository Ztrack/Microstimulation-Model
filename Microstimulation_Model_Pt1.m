clearvars; clc; close all;
load('InitialConditionsFull.mat')
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)

Stim_Center = 0; % If 1, uses a stimulus at center electrode instead of random
Num_Stim = 5;
I0 = 100;

%% Stimulation Locations

if Stim_Center == 1
    ElectrodeDist = sqrt((sx/2-ElectrodeX).^2 + (sy/2-ElectrodeY).^2);
    ElectrodeNo = find(ElectrodeDist == min(ElectrodeDist),1); % Finds the closest electrode to the center, stimulate only this electrode
else
    ElectrodeNo = randperm(100,Num_Stim);
end

Stim_Loc =  zeros(sx, sy);
for i = 1:length(ElectrodeNo)
Stim_Loc(ElectrodeY(ElectrodeNo(i))-Electroderadii:ElectrodeY(ElectrodeNo(i))+Electroderadii,ElectrodeX(ElectrodeNo(i))-Electroderadii:ElectrodeX(ElectrodeNo(i))+Electroderadii)  = 1;
end
Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
%% Extracellular Current Matrix


Ie_Axon_Neurons = sum(I0_Axon_Neurons(:,ElectrodeNo).*I0,2);
Ie_Soma_Neurons = sum(I0_Soma_Neurons(:,ElectrodeNo).*I0,2);
Ie_Soma_Axon_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons


Axon_Population_Matrix = Axon_Inhibitory_matrix + Axon_Excitatory_matrix; % Combined Axons from inhibitory + Excitatory neurons

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
    if Neuron_Type(i) == 1 % If it is an FS (Inhibitory)
        Lambda_Hat(i) = FS_Lambda + (FS_Lambda * Ie_Soma_Axon_Neurons(i)); % Lambda_Hat firing rate based off microstimulation
        Lambda(i) = FS_Lambda; % Lambda regular, no microstimulation
    end
end

Inhibitory_Effect = zeros(NumMotifs,1);
Inhibitory_Effect_Lambda = zeros(NumMotifs,1);
for i = 1:NumMotifs
    Inhibitory_Effect(i) = Lambda_Hat(Neuron_Inhibitory(i)).*Inhibitory_Factor; % Calculates the inhibitory effect from FS firing on RS neurons per motif
    Inhibitory_Effect_Lambda = FS_Lambda*Inhibitory_Factor;
end

for i = 1:NumNeurons % Applying inhibitory factor to firing rates of RS Neurons
    if Neuron_Type(i) == 2
        Lambda_Hat(i) = (RS_Lambda + (RS_Lambda * Ie_Soma_Axon_Neurons(i))) - Inhibitory_Effect(Neuron_Motif(i)); % Lambda hat firing rate based off stimulation & Inhibitory effect
        Lambda(i) = RS_Lambda - Inhibitory_Effect_Lambda; % Lambda regular, no microstimulation
    else
        % Inhib Neurons do not inhib themselves
    end
end

%% Poisson Spiking Model pt 1 - Lambda (No Stimulation)
Distances = zeros(length(ElectrodeNo),length(MotifX));
for i = 1:length(ElectrodeNo)
    Distances(i,:) = sqrt((ElectrodeX(ElectrodeNo(i))-MotifX).^2 + (ElectrodeY(ElectrodeNo(i)) - MotifY).^2); % Stores distances of every motif to every active electrode
end
Dist = sort(Distances(1,:));
Stimulated_Motifs = find(Distances(1,:) == Dist(1)); % Stores which motif is closest
Stimulated_Neurons = find(Neuron_Motif == Stimulated_Motifs(1));
    
NumTrials = 1000;
Oscillatory_Behavior = 1;
for i = 1:NumNeurons
    if intersect(Neuron_Inhibitory_Oscillatory,i) == i & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
        Lambda_Spikes = Oscillatory_PoissonGen(Lambda(i), dt, NumTrials);
    else % Otherwise use the simple function:
        Lambda_Spikes = Simple_PoissonGen(Lambda(i), dt, NumTrials);
    end
    Lambda_Hat_Spikes_Mean(i) =  mean(Lambda_Spikes);
end

Neuron_Pop_Spikes_Lambda = zeros(sx, sy);
for i = 1:NumNeurons
    Neuron_Pop_Spikes_Lambda(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii) = Lambda_Spikes(i);
end

%% Poisson Spiking Model pt 2 - Lambda Hat (Stimulated)
NumTrials = 1000;
Oscillatory_Behavior = 1;
for i = 1:NumNeurons
    if intersect(Neuron_Inhibitory_Oscillatory,i) == i & Oscillatory_Behavior == 1 % If the neuron has oscillatory behavior then use:
        Lambda_Hat_Spikes = Oscillatory_PoissonGen(Lambda_Hat(i), dt, NumTrials);
    else % Otherwise use the simple function:
        Lambda_Hat_Spikes = Simple_PoissonGen(Lambda_Hat(i), dt, NumTrials);
    end
    Lambda_Hat_Spikes_Mean(i) =  mean(Lambda_Hat_Spikes);
end


Neuron_Pop_Spikes = zeros(sx, sy);
for i = 1:length(NeuronX)
    Neuron_Pop_Spikes(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii) = Lambda_Hat_Spikes_Mean(i);
end

%% Plots for Model

figure;
map = [1 1 1; 0 0 0]; colormap(map);
imagesc(1./Stim_Distance); title('Active Electrode Locations'); xlabel('X Position (mm)'); ylabel('Y Position (mm)');

Neuron_Inhibitory_Population_Matrix = zeros(sx,sy);
Neuron_Excitatory_Population_Matrix = zeros(sx,sy);
for i = 1:NumNeurons
    if Neuron_Type(i) == 1
        Neuron_Inhibitory_Population_Matrix(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii) = 1;
    else
        if length(intersect(i,Neuron_Motion(1,:))) == 1 % If true, this neuron is a motion neuron
            Neuron_Excitatory_Population_Matrix(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii) = 3;
        else % Not a motion neuron
            Neuron_Excitatory_Population_Matrix(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii) = 2;
        end
    end
end
a = find(Neuron_Inhibitory_Population_Matrix == 1);
b = find(Neuron_Excitatory_Population_Matrix == 2);
c = find(Neuron_Excitatory_Population_Matrix == 3);
[y1,x1] = ind2sub([sx sy],a); [y2,x2] = ind2sub([sx sy],b); [y3,x3] = ind2sub([sx sy],c);
Axon_Population_Matrix(Axon_Population_Matrix > 0) = 1;

figure; 
map = [1 1 1];
imagesc(Axon_Population_Matrix); colormap(map);
hold on
a = find(Axon_Population_Matrix > 0);
[y4,x4] = ind2sub([sx sy],a);
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); %legend('Inhibitory','Excitatory');

figure; 
map = [1 1 1];
imagesc(Axon_Population_Matrix); colormap(map);
hold on
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Neural Population (Close Up)'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); %legend('Inhibitory','Excitatory');
range = [2000 3000]; % Zoomed in range for plot below
xlim(range); ylim(range);

PadMatrix = zeros(sx,sy); px=sx/5; py=sy/3;
for i = 1:NumPads
    rangex = ((Padloc(i,3)-1)*px)+1:(Padloc(i,3)*px)-1;
    rangey = ((Padloc(i,2)-1)*py)+1:(Padloc(i,2)*py)-1;
    PadMatrix(rangey,rangex) = i;
end
map = [0 0 0; 239 122 37; 81 182 106; 50 127 183; 152 85 159; 225 152 153;]./255; % Colormap. First number is 0, then 1 to 15.
c1 = .85;
c2 = .7;
map = [map ; map(2,1:3)*c1 ; map(3,1:3)*c1 ; map(4,1:3)*c1 ; map(5,1:3)*c1 ; map(6,1:3)*c1];
map = [map ; 1 1 1 ; map(3,1:3)*c2 ; map(4,1:3)*c2 ; map(5,1:3)*c2 ; map(6,1:3)*c2];

figure; imagesc(PadMatrix); colormap(map); title('Pads'); % Color coded map representing different pads

figure; imagesc(PadMatrix); colormap(map); title('Pads with electrodes'); % Color coded map representing different pads
hold on;
plot(ElectrodeX,ElectrodeY,'.','color','black');


figure; imagesc(PadMatrix); colormap(map); title('Neuron Synaptic Connections'); % Color coded map representing different pads
hold on
for i = 1:length(Neuron_Connected)
        x1 = NeuronX(Neuron_Connected(i,1)); x2 = NeuronX(Neuron_Connected(i,2));
        y1 = NeuronY(Neuron_Connected(i,1)); y2 = NeuronY(Neuron_Connected(i,2));
        line([x1,x2], [y1,y2],'Color', 'black');
end

%%

% Stimulated motifs / neuron pop maps:
rangex = Rangex(Stimulated_Motifs,:); rangey = Rangey(Stimulated_Motifs,:);
figure; imagesc(Neuron_Pop_Spikes(rangey,rangex)); title('Stimulated Motif'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 max(Lambda_Hat(Stimulated_Neurons))]);
hold on
a = find(Axon_Population_Matrix(rangey,rangex) > 0);
[y1,x1] = ind2sub([length(rangex) length(rangey)],a);
plot(x1,y1,'.','color','Black');
grad=colorGradient([0 0 1],[1 1 0],128);
map = [1 1 1; grad];
colormap(map);



figure; imagesc(Neuron_Pop_Spikes); title('Neuron Response Stimulated'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 100]);
colormap(map);

% Non-Stimulated motifs / neuron pop maps:
figure; imagesc(Neuron_Pop_Spikes_Lambda(rangey,rangex)); title('Non-Stimulated Motif'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 max(Lambda_Hat(Stimulated_Neurons))]);
hold on
a = find(Axon_Population_Matrix(rangey,rangex) > 0);
[y1,x1] = ind2sub([length(rangex) length(rangey)],a);
plot(x1,y1,'.','color','Black');
colormap(map);

figure; imagesc(Neuron_Pop_Spikes_Lambda); title('Neuron Response Non-Stimulated'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar;
colormap(map);

% motifs difference map:
figure; imagesc(Neuron_Pop_Spikes(rangey,rangex) - Neuron_Pop_Spikes_Lambda(rangey,rangex)); title('Motif Difference Map'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([0 max(Lambda_Hat(Stimulated_Neurons)-max(Lambda(Stimulated_Neurons)))]);
hold on
a = find(Axon_Population_Matrix(rangey,rangex) > 0);
[y1,x1] = ind2sub([length(rangex) length(rangey)],a);
plot(x1,y1,'.','color','Black')
colormap(map);

figure; imagesc(Neuron_Pop_Spikes-Neuron_Pop_Spikes_Lambda); title('Neuron Response Difference'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); colorbar; caxis([1 100]);
colormap(map);

%% Functions
function [MotifNeuronX,MotifNeuronY,Axon_Loc_Inhibitory,Neuron_Loc] = MotifCreator(Mx,My,NumNeuronsMotif,NeuronRadii,Neuron_spacing)

for i = 1:NumNeuronsMotif
    if i == 1 % If i is an inhibitory neuron
        MotifNeuronX(i) = randi([NeuronRadii+1,Mx-NeuronRadii-1],1); % x coordinates
        MotifNeuronY(i) = NeuronRadii+1; % y coordinates
    else % if i is an Excitatory neuron, for i = 2:5
        Counter = 2;
        while Counter~=NumNeuronsMotif+1
            TestX = randi([NeuronRadii+1,Mx-NeuronRadii-1],1); % X value we are testing
            TestY = randi([MotifNeuronY(1)+Neuron_spacing*2,Mx-NeuronRadii-1],1); % Y value we are testing
            Distances = sqrt((TestX-MotifNeuronX).^2 + (TestY - MotifNeuronY).^2); % Stores distances to every other neuron
            MinDistance = min(Distances); % Finds the closest neuron to this point
            XDististances = sqrt((TestX-MotifNeuronX).^2);
            MinXDistance = min(XDististances);
            if MinDistance > Neuron_spacing && MinXDistance > Mx/10 % If the distance far enough away, it is accepted
                MotifNeuronX(Counter) = TestX;
                MotifNeuronY(Counter) = TestY;
                Counter = Counter + 1;
            end
        end
    end
end
Neuron_Loc_point =  zeros(Mx, My);
Neuron_Loc = zeros(Mx, My);
for i = 1:length(MotifNeuronX)
    Neuron_Loc_point(MotifNeuronY(i),MotifNeuronX(i)) =  1;
end

% Axon Location Matrix
Axon_Loc_Inhibitory =  zeros(Mx, My);
for i = 1:NumNeuronsMotif+1
    if i == 1 % If i is an inhibitory neuron
        xc = [MotifNeuronX(i) MotifNeuronX(i)]; % x coordinates
        yc = [MotifNeuronY(i) MotifNeuronY(i)+Neuron_spacing]; % y coordinates
    elseif i >1 && i<=NumNeuronsMotif % here we draw the axons going up to a Y coordinate just below the inhib neuron
        xc = [MotifNeuronX(i) MotifNeuronX(i)]; 
        yc = [MotifNeuronY(i) MotifNeuronY(1)+Neuron_spacing];
    elseif i == NumNeuronsMotif+1 % Here we draw a line to connect the neuron motif
        xc = [min(MotifNeuronX) max(MotifNeuronX)];
        yc = [MotifNeuronY(1)+Neuron_spacing MotifNeuronY(1)+Neuron_spacing];
    end
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind(size(Axon_Loc_Inhibitory), rIndex, cIndex); % Linear indices
    Axon_Loc_Inhibitory(lindex) = 1;  % Gives us all points where the axons cross on x,y plane
end

for i = 1:Mx     
    for j = 1:My
        if Neuron_Loc_point(i,j) == 1
                Axon_Loc_Inhibitory(i-NeuronRadii:i+NeuronRadii, j-NeuronRadii:j+NeuronRadii) = 0; % Refining Axon_Loc_Inhibitory to not include spaces occupied by neurons
                Neuron_Loc(i-NeuronRadii:i+NeuronRadii, j-NeuronRadii:j+NeuronRadii) = 1; % All spaces occupied by neurons in motif
        end
    end
end

end

function [rp_x,rp_y] = XYFILTER(rp_x,rp_y,Mx,My,PadsXLength,PadsYLength,ii)
% First we filter the X values
if ii == 1 || ii == 6 || ii == 11 % Left Column
    for i = 1:length(rp_x)
        if rp_x(i) < Mx
            rp_x(i) = NaN;
            rp_y(i)= NaN;
        end
    end
elseif ii == 5 || ii == 10 || ii == 15 % Right Column
    for i = 1:length(rp_x) 
        if rp_x(i) > PadsXLength-Mx
            rp_x(i) = NaN;
            rp_y(i)= NaN;
        end
    end
else
    % If middle row, do nothing
end

% Second we filter the Y values
if ii == 1 || ii == 2 || ii == 3 || ii == 4 || ii == 5 % Top Row
    for i = 1:length(rp_y)
        if rp_y(i) < My
            rp_x(i) = NaN;
            rp_y(i)= NaN;
        end
    end
elseif ii == 11 || ii == 12 || ii == 13 || ii == 14 || ii == 15 % Bottom row
    for i = 1:length(rp_y)
        if rp_y(i) > PadsYLength-My
            rp_x(i) = NaN;
            rp_y(i)= NaN;
        end
    end
else
    % If middle row, do nothing
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