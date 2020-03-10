clear all
%%
Res = 1; % Resolution factor. 1 = 1um/square. 2 = .5 um/square, etc.
X_length = 4800; % Size of work area. Must be divisible by 3 & 5 (Number of pads in X,Y)
Y_length = 4800;
NeuronRadii = 5; % in micrometers
Neuron_spacing = 40; % micrometers seperation of neurons. Too large a value can cause errors
Num_Stim = 5; % Number of stimulations in present in matrix area
T = 1; % Duration of experiment in seconds
I0 = 100; % Amps of current injected
sx = X_length*Res;
sy = Y_length*Res;
Mx = 200; % Motif X length
My = 200; % Motif Y Length
NumNeuronsMotif = 10; % # Neuron in Motif. Must be multiples of 5 to satisfy ratio. Every 1st neuron is Inhibitory
NumMotifs = 120; % Number of motifs
NumNeurons = NumMotifs*NumNeuronsMotif;
NumPadsX = 5; % Number of pads to divide X_Length
NumPadsY = 3;
NumPads = NumPadsX*NumPadsY; % Total number of pads in area
NumMotifsPerPad = NumMotifs/NumPads; % Number of motifs per pad
PadsXLength = X_length/NumPadsX; % Length of a X pad
PadsYLength = Y_length/NumPadsY; % Length of a Y pad
Inhibitory_Factor = .01; % 0 to 1 Greater factor will inhibit RS neurons more

%% Neuron Definitions

Neuron_Type = repmat([1,1,1,1,1,2,2,2,2,2],1,NumMotifs); % Defines what type of neuron this number is. Inhibitory =1, Excitatory = 2
Neuron_Motif = repmat(1:1:NumMotifs,NumNeuronsMotif,1); Neuron_Motif = Neuron_Motif(:)'; % Defines which motif the neuron belongs to

Neuron_Inhibitory = find(Neuron_Type == 1); % Array with all inhibitory neurons
Neuron_Excitatory = find(Neuron_Type == 2); % Array with all Excitatory neurons
Pad_Motif = repmat(1:NumPads,NumMotifs/NumPads,1); Pad_Motif = Pad_Motif(:); % Stores which pad the motif belongs to
Pad_Neuron = repmat(1:NumPads,NumNeurons/NumPads,1); Pad_Neuron = Pad_Neuron(:); % Stores which pad the motif belongs to
Neuron_Motion = sort(Neuron_Excitatory(randperm(length(Neuron_Excitatory),NumNeurons*0.3))); % Neurons selected for motion
Neuron_NonMotion = ones(NumNeurons,1); Neuron_NonMotion(Neuron_Motion) = 0; Neuron_NonMotion = find(Neuron_NonMotion == 1); % Non-Motion neurons
Neuron_Motion(2,:) = randi([1,5],[1,length(Neuron_Motion)]); % Gives every motion neuron an orientation
Coding_Probability = [6/13 3/16 4/22 4/16 2/5]; % Connection Probability, given by Ko et al 2011
Neuron_Connected = []; k =1;
for i = 1:5
    Coding = Neuron_Motion(:,2) == i; % Finds indices where neurons are coded for same motion
    Coding = Neuron_Motion(Coding,1); % Converts to Neuron Number
    for j = 1:length(Coding)
        X_random = rand(length(Coding),1);
        if Coding_Probability(i) > X_random(j) % If the probability of connection is greater we store it
            Neuron_Connected(k,1) = Coding(j);
            Neuron_Connected(k,2) = Coding(randi([1 length(Coding)],[1 1])); % Connects to a random neuron in same coding
            k = k+1;
        end
    end
end
Neuron_Inhibitory_Oscillatory = Neuron_Inhibitory(randperm(numel(Neuron_Inhibitory), ceil(length(Neuron_Inhibitory)*.33)));

%% Motif Location Selector
Padloc = [1,1,1;2,1,2;3,1,3;4,1,4;5,1,5;6,2,1;7,2,2;8,2,3;9,2,4;10,2,5;11,3,1;12,3,2;13,3,3;14,3,4;15,3,5;]; % #pad, column, row - Location of pad
Neuron_Population_Matrix = zeros(sx, sy); % Initialize stimulation location matrix size nxm
MotifX = zeros(NumMotifs,1); MotifY = zeros(NumMotifs,1); % Defines center point of motif in the larger population matrix
PadSize = zeros(PadsYLength, PadsXLength);
kk = 1;
for ii = 1:NumPads
    rp = randperm(numel(PadSize)); % Randomly orders the matrix in indices form
    [rp_x,rp_y] = ind2sub([PadsXLength PadsYLength], rp);% Returns the x,y of the random permutation
    [rp_x,rp_y] = XYFILTER(rp_x,rp_y,Mx,My,PadsXLength,PadsYLength,ii); % Filter so motif does not crop out of defined area % To make this modular, add NumPadMax,NumPadMin for x,y and ,a
    rp_x = rp_x + PadsXLength*(Padloc(ii,3)-1);
    rp_y = rp_y + PadsYLength*(Padloc(ii,2)-1);
    Counter = 1;
    k = 1;
    while Counter~=NumMotifsPerPad+1 % While loop tests distance to other motif locations
        TestX = rp_x(k); % X value we are testing
        TestY = rp_y(k); % Y value we are testing
        Distances = sqrt((TestX-MotifX).^2 + (TestY - MotifY).^2); % Stores distances to every other motif
        MinDistance = min(Distances); % Finds the closest neuron to this point
        if MinDistance > Mx*sqrt(2) % If the distance far enough away, it is accepted
            MotifX(kk) = TestX;
            MotifY(kk) = TestY;
            Counter = Counter + 1;
            k = k + 1;
            kk = kk + 1;
        else
            k= k + 1;
        end
    end
end

%% Neural Population Matrix

Axon_Inhibitory_matrix = zeros(sx,sy);
Neuron_Points_Matrix = zeros(sx,sy);
Neuron_Indices = zeros(NumNeurons,1);
NeuronX = zeros(NumNeurons,1); NeuronY = zeros(NumNeurons,1); % X,Y locations of each neuron
Rangex = zeros(NumMotifs,Mx); Rangey = zeros(NumMotifs,My); % Range of x,y values in which our motifs exist
Orientation = randi([1,4],NumMotifs,1); % Random number 1-4 for Orientation of motif. 1 = 90, 2 = 180, 3 = 270, 4 = 360
kk = 1; % # Neuron #
for i = 1:NumMotifs
    [MotifNeuronX,MotifNeuronY,Axon_Loc_Inhibitory,Neuron_Loc] = MotifCreator(Mx,My,NumNeuronsMotif,NeuronRadii,Neuron_spacing); % Motif Generator function
    Rangex(i,:) = MotifX(i)-length(Neuron_Loc)/2:MotifX(i)+length(Neuron_Loc)/2-1; % Location of motif start to end x
    Rangey(i,:) = MotifY(i)-length(Neuron_Loc)/2:MotifY(i)+length(Neuron_Loc)/2-1; % Location of motif start to end y
    for j = 1:length(MotifNeuronX) % Rotating NeuronX,NeuronY
        Neuron_points_motif = zeros(size(Neuron_Loc));
        Neuron_points_motif(MotifNeuronY(j),MotifNeuronX(j)) = 1; % Rotating one neuron at a time to accurately determine it's location
        for k = 1:Orientation(i)
            Neuron_points_motif = rot90(Neuron_points_motif);
        end
        [Y,X] = ind2sub(size(Neuron_points_motif),find(Neuron_points_motif ==1)); % New values for Motif Neuron location
        NeuronX(kk) = min(Rangex(i,:)) + X; NeuronY(kk) = min(Rangey(i,:)) + Y; % Translating new motif neuron location to larger matrix location
        kk = kk + 1; % For every neuron
    end
    for j = 1:Orientation(i) % Rotates motif neuron locations and axon locations
        Axon_Loc_Inhibitory = rot90(Axon_Loc_Inhibitory); % Rotating inhibitory axon network
        Neuron_Loc = rot90(Neuron_Loc); % Rotating neuron map
    end
    Neuron_Population_Matrix(Rangey(i,:),Rangex(i,:))=Neuron_Loc; % Constructs neuron pop matrix from each motif
    Axon_Inhibitory_matrix(Rangey(i,:),Rangex(i,:))=Axon_Loc_Inhibitory; % Constructs inhibitory axon pop from each motif
end
for i = 1: NumNeurons
    Neuron_Points_Matrix(NeuronY(i),NeuronX(i))= 1; % Builds neuron points matrix, where a neuron is defined by their center point
end

%% Stimulation Locations Matrix

Electroderadii = 10;
ElectrodeX = []; ElectrodeY = [];
for i = 1:10 % Number of electrodes = 100, 10x10 over 4mm in center of map. Does not create an electrode at percent center!
    for j = 1:10
        ElectrodeX(i,j) = (i-1) * sx/10 + sy/20;
        ElectrodeY(i,j) = (j-1) * sy/10 + sy/20;
    end
end
ElectrodeX = ElectrodeX(:); ElectrodeY = ElectrodeY(:);

Stim_Distance_Map = zeros(length(ElectrodeX),sx,sy);
for i = 1:length(ElectrodeX)
    Stim_Loc =  zeros(sx, sy); % Location of stimulus
    Stim_Loc(ElectrodeY(i)-Electroderadii:ElectrodeY(i)+Electroderadii,ElectrodeX(i)-Electroderadii:ElectrodeX(i)+Electroderadii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
    Stim_Distance_Map(i,:,:) = Stim_Distance;
end

%% Optogenetics

Optogenetics = 0;

if Optogenetics == 1
    Opto_Electrodeno = 1; % Number of electrode emitting light
    Optop_Loc(ElectrodeY(Opto_Electrodeno)-Electroderadii:ElectrodeY(Opto_Electrodeno)+Electroderadii,ElectrodeX(Opto_Electrodeno)-Electroderadii:ElectrodeX(Opto_Electrodeno)+Electroderadii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Opto_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
end

%% Neuron And Axon Area Summation

Axon_Excitatory_matrix = zeros(sx,sy);

Excitatory_Axon_Storage = nan(NumNeurons,500); % Stores indices of this neuron's axon
I0_Soma_Neurons = zeros(NumNeurons,length(ElectrodeX));
I0_Axon_Neurons = zeros(NumNeurons,length(ElectrodeX)); % Stores 1/r^2 area component of Axons

for i = 1:NumNeurons
    I0_Axon_NeuronsA = zeros(1,length(ElectrodeX));
    I0_Soma_NeuronsA = zeros(1,length(ElectrodeX));
    if Neuron_Type(i) == 1
        x = MotifX(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifX(Neuron_Motif(i))+length(Neuron_Loc)/2-1; % Location of motif does not change with Orientation
        y = MotifY(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifY(Neuron_Motif(i))+length(Neuron_Loc)/2-1;
        for j = 1:length(ElectrodeX)
                Stim_Distance = squeeze(Stim_Distance_Map(j,:,:));
                I0_Axon_NeuronsA(j) = sum(sum(Axon_Inhibitory_matrix(y,x)./(Stim_Distance(y,x).^2))); % sums all inhibitory axon current feedback
                I0_Soma_NeuronsA(j) = sum(sum(1./Stim_Distance(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma
        end
    elseif Neuron_Type(i) == 2 % Creating a new axon matrix for Excitatory neurons. Must be made one at a time to calculate current effect on every neuron
        Axon_Excitatory_matrix_single = zeros(sx,sy); % Axon matrix for a single neuron
            re = [0,0,0,0]; re(randi([1,4],1)) = 1; % Random element selection to choose orientation of axon
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
                I0_Axon_NeuronsA(j) = sum(sum(1./(Stim_Distance(lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
                I0_Soma_NeuronsA(j) = sum(sum(1./Stim_Distance(NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma
            end
            Axon_Excitatory_matrix = Axon_Excitatory_matrix + Axon_Excitatory_matrix_single; % builds axon excitatory matrix
    end
    I0_Axon_Neurons(i,:) = I0_Axon_NeuronsA;
    I0_Soma_Neurons(i,:) = I0_Soma_NeuronsA;
end
    


clear Stim_Distance_Map;
save InitialConditionsFullB.mat
%% Functions
function [MotifNeuronX,MotifNeuronY,Axon_Loc_Inhibitory,Neuron_Loc] = MotifCreator(Mx,My,NumNeuronsMotif,NeuronRadii,Neuron_spacing)

for i = 1:NumNeuronsMotif
    if i == 1 || 2 || 3 || 4 || 5 % If i is an inhibitory neuron
        MotifNeuronX(i) = 25+15*i; % x coordinates
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
    if i <= 5 % If i is an inhibitory neuron
        xc = [MotifNeuronX(i) MotifNeuronX(i)]; % x coordinates
        yc = [MotifNeuronY(i) MotifNeuronY(i)+Neuron_spacing]; % y coordinates
    elseif i <= 10  && i<=NumNeuronsMotif % here we draw the axons going up to a Y coordinate just below the inhib neuron
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