clear all;
%%
Res = 1; % Resolution factor. 1 = 1um/square. 2 = .5 um/square, etc.
X_length = 4800; % Size of work area. Must be divisible by 3 & 5 (Number of pads in X,Y)
Y_length = 4800;
NeuronRadii = 5; % in micrometers
Neuron_spacing = 40; % micrometers seperation of neurons. Too large a value can cause errors
Num_Stim = 5; % Number of stimulations in present in matrix area
T = 1; % Duration of experiment in seconds
dt = T/1000;
I0 = 100; % Amps of current injected
sx = X_length*Res;
sy = Y_length*Res;
Mx = 200; % Motif X length
My = 200; % Motif Y Length
NumNeuronsMotif = 5; % # Neuron in Motif. Must be multiples of 5 to satisfy ratio. Every 1st neuron is Inhibitory
NumMotifs = 120; % Number of motifs
NumNeurons = NumMotifs*NumNeuronsMotif;
NumPadsX = 5; % Number of pads to divide X_Length
NumPadsY = 3;
NumPads = NumPadsX*NumPadsY; % Total number of pads in area
NumMotifsPerPad = NumMotifs/NumPads; % Number of motifs per pad
PadsXLength = X_length/NumPadsX; % Length of a X pad
PadsYLength = Y_length/NumPadsY; % Length of a Y pad
Inhibitory_Factor = .01; % 0 to 1 Greater factor will inhibit RS neurons more
Inhibitory_Oscillatory_Percentage = .33; % Percentage of inhibitory neurons which oscillate

%% Neuron Definitions

Neuron_Type = repmat([1,2,2,2,2],1,NumMotifs); % Defines what type of neuron this number is. Inhibitory =1, Excitatory = 2
Neuron_Motif = repmat(1:1:NumMotifs,NumNeuronsMotif,1); Neuron_Motif = Neuron_Motif(:)'; % Defines which motif the neuron belongs to

Neuron_Inhibitory = find(Neuron_Type == 1); % Array with all inhibitory neurons
Neuron_Excitatory = find(Neuron_Type == 2); % Array with all Excitatory neurons
Pad_Motif = repmat(1:NumPads,NumMotifs/NumPads,1); Pad_Motif = Pad_Motif(:); % Stores which pad the motif belongs to
Pad_Neuron = repmat(1:NumPads,NumNeurons/NumPads,1); Pad_Neuron = Pad_Neuron(:); % Stores which pad the motif belongs to
Neuron_Motion = sort(Neuron_Excitatory(randperm(length(Neuron_Excitatory),NumNeurons*0.3))); % Neurons selected for motion
Neuron_NonMotion = ones(NumNeurons,1); Neuron_NonMotion(Neuron_Motion) = 0; Neuron_NonMotion = find(Neuron_NonMotion == 1); % Non-Motion neurons
Neuron_Motion(2,:) = randi([1,2],[1,length(Neuron_Motion)]); % Gives every motion neuron an orientation 0 or 90 degrees
Coding_Probability = [10/26 4/24]; % Connection Probability, given by Ko et al 2011
Directional_Probability = [.1 .1; .4 .1]; % Unidirectional 0 degrees = 10%, 90 degrees = 10% ; Bidirectional 0 degrees = 40%, 90 degrees = 10%
Neuron_Connected = []; k =1;
for i = 1:length(Neuron_Motion)
    for j = 1:length(Neuron_Motion)
        if Coding_Probability(Neuron_Motion(2,i)) > rand(1) && i~=j
                  Neuron_Connected(k,1) = Neuron_Motion(1,i); % Source Neuron
                  Neuron_Connected(k,2) = Neuron_Motion(1,j); % Connected Neuron
                  if Coding_Probability(Neuron_Motion(2,i)) == 10/26
                      Neuron_Connected(k,3) = (.8 > rand(1)); % 1 = Bidirectional, 0 = Unidirectional in direction of connected pair
                  elseif Coding_Probability(Neuron_Motion(2,i)) == 4/24
                      Neuron_Connected(k,3) = (.5 > rand(1));
                  end
                  k = k+1;
        end
    end
end

Neuron_Inhibitory_Oscillatory = Neuron_Inhibitory(randperm(numel(Neuron_Inhibitory), ceil(length(Neuron_Inhibitory)*Inhibitory_Oscillatory_Percentage)));
Neuron_Inhibitory_Non_Oscillatory = setdiff(Neuron_Inhibitory,Neuron_Inhibitory_Oscillatory);
Neuron_Type2 = ones(1,600); Neuron_Type2(Neuron_Inhibitory_Non_Oscillatory) = 2; Neuron_Type2(Neuron_Inhibitory_Oscillatory) = 3;

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

%% Neuron And Axon Area Summation

Axon_Excitatory_matrix = zeros(sx,sy);

Excitatory_Axon_Storage = nan(NumNeurons,500); % Stores indices of this neuron's axon
I0_Soma_Neurons = zeros(NumNeurons,length(ElectrodeX));
I0_Axon_Neurons = zeros(NumNeurons,length(ElectrodeX)); % Stores 1/r^2 area component of Axons
Axon_Direction = randi([1,4],NumNeurons,1); % 1 = x value down (left),2 = Y value down (up), 3 = X value up (right), 4 = Y value up (down)
for i = 1:length(Neuron_Inhibitory)
    Axon_Direction(Neuron_Inhibitory(i)) = Orientation(i); % Sets stored value for inhibitory axon direction
end
Axon_Indices = cell(NumNeurons,1);
Soma_Indices = cell(NumNeurons,1);
for i = 1:NumNeurons
Soma_Indices(i) = {sub2ind([sx sy],NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii)};
end

for i = 1:NumNeurons
    I0_Axon_Neurons1 = zeros(1,length(ElectrodeX));
    I0_Soma_Neurons1 = zeros(1,length(ElectrodeX));
    if Neuron_Type(i) == 1
        x = MotifX(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifX(Neuron_Motif(i))+length(Neuron_Loc)/2-1; % Location of motif does not change with Orientation
        y = MotifY(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifY(Neuron_Motif(i))+length(Neuron_Loc)/2-1;
        lindex = find(Axon_Inhibitory_matrix(y,x)==1);
        for j = 1:length(ElectrodeX)
                I0_Axon_Neurons1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % sums all inhibitory axon current feedback
                I0_Soma_Neurons1(j) = sum(sum(1./Stim_Distance_Map(j,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma % Summation 1/r^2 area component of soma
        end
        
    elseif Neuron_Type(i) == 2 % Creating a new axon matrix for Excitatory neurons. Must be made one at a time to calculate current effect on every neuron
            re = [0,0,0,0]; re(Axon_Direction(i)) = 1; % Random element selection to choose orientation of axon
            rl = randi([50+NeuronRadii,480+NeuronRadii],1); % Random length of axon + NeuronRadii
            xc = [NeuronX(i)+(re(3)*NeuronRadii)-(re(1)*NeuronRadii) NeuronX(i)+(re(3)*rl)-(re(1)*rl)]; % x coordinates
            yc = [NeuronY(i)+(re(4)*NeuronRadii)-(re(2)*NeuronRadii) NeuronY(i)+(re(4)*rl)-(re(2)*rl)]; % y coordinates
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
            lindex = sub2ind([sx sy], rIndex,cIndex); % Linear indices
            for j = 1:length(ElectrodeX)
                I0_Axon_Neurons1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
                I0_Soma_Neurons1(j) = sum(sum(1./Stim_Distance_Map(j,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma
            end
            
    end
    Axon_Indices{i,:} = lindex(:);
    I0_Axon_Neurons(i,:) = I0_Axon_Neurons1;
    I0_Soma_Neurons(i,:) = I0_Soma_Neurons1;
end

% Axon Generation for motion-tuned neural connections
I0_Motion_Neuron_Connections = zeros(length(Neuron_Connected),length(ElectrodeX));
for i = 1:length(Neuron_Connected)
    I0_Motion_Neuron_Connections1 = zeros(1,length(ElectrodeX));
    xc = [NeuronX(Neuron_Connected(i,1)) NeuronX(Neuron_Connected(i,2))]; % x coordinates
    yc = [NeuronY(Neuron_Connected(i,1)) NeuronY(Neuron_Connected(i,2))]; % y coordinates
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind([sx sy], rIndex,cIndex); % Linear indices
    for j = 1:length(ElectrodeX)
        I0_Motion_Neuron_Connections1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
    end
    I0_Motion_Neuron_Connections(i,:) = I0_Motion_Neuron_Connections1(j);
end
I0_Motion_Neurons = zeros(NumNeurons,length(ElectrodeX));
for i = 1:length(I0_Motion_Neuron_Connections)
    I0_Motion_Neurons(Neuron_Connected(i,1),:) = I0_Motion_Neurons(Neuron_Connected(i,1),:) + I0_Motion_Neuron_Connections(i,:);
    I0_Motion_Neurons(Neuron_Connected(i,2),:) = I0_Motion_Neurons(Neuron_Connected(i,2),:) + I0_Motion_Neuron_Connections(i,:);
end

for i = 1:NumNeurons
Axon_Excitatory_matrix(Axon_Indices{i,:}) = 1; % Builds excitatory axon matrix
end

Directional_Current_Mult = zeros(NumNeurons,length(ElectrodeX));
Axon_Direction_Theta = zeros(NumNeurons,1);
if Axon_Direction(i) == 1 %(left = 180 degrees)
    Axon_Direction_Theta(i) = 180;
elseif Axon_Direction(i) == 2 % (up = 90 degrees)
    Axon_Direction_Theta(i) = 90;
elseif Axon_Direction(i) == 3 % (right = 0 degrees)
    Axon_Direction_Theta(i) = 0;
elseif Axon_Direction(i) == 4 % (down, = -90 or 270 degrees)
    Axon_Direction_Theta(i) = -90;
end
theta_threshold = 45; % angle difference maximum

for i = 1:NumNeurons
    x2 = NeuronX(i);
    y2 = NeuronY(i);
    
for ii = 1:length(ElectrodeX)
    
    x1 = ElectrodeX(ii);
    y1 = ElectrodeY(ii);
    theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
    normDeg = mod(Axon_Direction_Theta(i)-theta,360);
    absDiffDeg = min(360-normDeg, normDeg);
    
    if absDiffDeg < theta_threshold % If angle meets critera then:
        Directional_Current_Mult(i,ii) = 1; % Stores directional current multiplier as full
    else % If angle is not within criteria
        Directional_Current_Mult(i,ii) = 0.25; % Stores current summation as depreciated
    end
end
end

%% Cleaning output

clear ('Stim_Distance_Map','rp','rp_x','rp_y','Stim_Distance','Stim_Loc','PadSize','Neuron_points_Matrix','Ed');
save InitialConditionsFull.mat;
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