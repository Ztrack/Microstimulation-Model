clearvars; clc;
%% IC

X_length = 4800; % Size of work area. Must be divisible by 3 & 5 (Number of pads in X,Y)
Y_length = 4800;
NeuronRadii = 5; % in micrometers
T = 1; % Duration of experiment in seconds
sx = X_length;
sy = Y_length;


NumNeurons = 980;
NeuronY(1:NumNeurons) = sy/2; NeuronX = [];
j = 1;
for i = 1:NumNeurons/2
    NeuronX(i) = sx/2 + 4*j;
    j = j+1;
end
NeuronX = repmat(NeuronX,1,2);
Neuron_Type = ones(size(NeuronX))+1;
Neuron_Pair = 1:NumNeurons/2; Neuron_Pair = repmat(Neuron_Pair,1,2);
Axon_Orientation(1:NumNeurons/2) = 1; Axon_Orientation((NumNeurons/2)+1:NumNeurons) = 2;
Axon_Toward = find(Axon_Orientation == 1); % Axons going toward field
Axon_Away = find(Axon_Orientation == 2); % Axons going away from electrical field

%% Stimulation Locations Matrix

Electroderadii = 10;
Stim_Loc =  zeros(sx, sy); % Location of stimulus
ElectrodeX = sx/2; ElectrodeY = sy/2;
Stim_Loc(ElectrodeY-Electroderadii:ElectrodeY+Electroderadii,ElectrodeX-Electroderadii:ElectrodeX+Electroderadii) = 1; % Builds electrode stimulus location matrix
Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
Stim_Distance_Map(1,:,:) = Stim_Distance;

%% Neuron And Axon Area Summation

I0_Soma_Neurons = zeros(NumNeurons,length(ElectrodeX));
I0_Axon_Neurons = zeros(NumNeurons,length(ElectrodeX)); % Stores 1/r^2 area component of Axons

for i = 1:NumNeurons
    Axon_Excitatory_matrix_single = zeros(sx,sy); % Axon matrix for a single neuron
    if Axon_Orientation(i) == 1
        re = [0,1,0,0]; % Toward Field
    else
        re = [1,0,0,0]; % Away from field
    end
    rl = randi([50+NeuronRadii,480+NeuronRadii],1); % Random length of axon + NeuronRadii
    xc = [NeuronX(i)+(re(1)*NeuronRadii)-(re(2)*NeuronRadii) NeuronX(i)+(re(1)*rl)-(re(2)*rl)]; % x coordinates
    yc = [NeuronY(i)+(re(3)*NeuronRadii)-(re(4)*NeuronRadii) NeuronY(i)+(re(3)*rl)-(re(4)*rl)]; % y coordinates
    if xc(2) <= 0 % Simple filter Makes sure the axon end does not go out of bounds in x
        xc(2) = 1;
    elseif xc(2) >= sx
        xc(2) = sx;
    end
    if yc(2) <= 0 %Simple filter Makes sure the axon end does not go out of bounds in y
        yc(2) = 1;
    elseif yc(2) >= sy
        yc(2) = sy;
    end
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind(size(Axon_Excitatory_matrix_single), rIndex,cIndex); % Linear indices
    for j = 1:length(ElectrodeX)
        Stim_Distance = squeeze(Stim_Distance_Map(j,:,:));
        I0_Axon_Neurons(i,j) = sum(sum(1./(Stim_Distance(lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
        I0_Soma_Neurons(i,j) = sum(sum(1./Stim_Distance(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma
    end
end

save InitialConditions_Pt2.mat