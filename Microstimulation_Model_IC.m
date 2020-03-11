clear all;

%% Initial Condition Parameters

% Neuron properties
neuron.radii = 5; % in micrometers
Electrode.radii = 10; % in micrometers
Inhibitory_Factor = .01; % 0 to 1 Greater factor will inhibit RS neurons more
theta_threshold = 45; % angle difference threshold - If the neuron axon is out of phase by at least this much to the current-stimulus, the electric field acting on the neuron is weakened to 25%.

% Population Properties
NumNeurons = 1000; % Must be multiples of 5 to satisfy ratio, if using 1:4 ratio. Every 1st neuron is Inhibitory
Excitatory_Inhibitory_Ratio = 1/4; % % Number of Inhibitory / Excitatory Neurons in a motif
NumNeuronsPerMotif = 5;
NumMotifs = NumNeurons/NumNeuronsPerMotif; % # Neuron in Motif.
Inhibitory_Oscillatory_Percentage = .33; % Percentage of inhibitory neurons which oscillate
NeuronMotionRatio = 0.2; % Ratio of Motion Neurons to non-motion neurons. based on population of excitatory neurons

% Spatial Properties
sx = 4800; % Size of work area
sy = 4800; % Size of work area
NumPadsX = 5; % Number of pads in X Direction
NumPadsY = 3; % Number of pads in Y direction
PadsXLength = sx/NumPadsX; % Length of a X pad
PadsYLength = sy/NumPadsY; % Length of a Y pad
NumPads = NumPadsX*NumPadsY; % Total number of pads in sx*sy area
PadSize = zeros(PadsYLength, PadsXLength);

% Temporal Properties
T = 1; % Duration of experiment in seconds
dt = T/1000; % Time step length

%% Neuron Definitions

neuron.number = 1:NumNeurons;
neuron.type = repmat([1,2,2,2,2],1,NumMotifs); % Defines what type of neuron this number is. Inhibitory =1, Excitatory = 2
neuron.motif = repmat(1:1:NumMotifs,NumNeuronsMotif,1); neuron.motif = neuron.motif(:)'; % Defines which motif the neuron belongs to
neuron.inhibitory = find(neuron.type == 1); % Array with all inhibitory neurons
neuron.excitatory = find(neuron.type == 2); % Array with all Excitatory neurons
neuron.motion.number = sort(neuron.excitatory(randperm(length(neuron.excitatory),floor(NumNeurons*NeuronMotionRatio)))); % Neurons selected for motion
neuron.nonMotion.number = ones(NumNeurons,1); neuron.nonMotion.number(neuron.motion.number) = 0; neuron.nonMotion.number = find(neuron.nonMotion.number == 1); % Non-Motion neurons
neuron.motion.direction = randi([1,2],[1,length(neuron.motion.number)]); % Gives every motion neuron an orientation 0 or 90 degrees
Coding_Probability = [10/26 4/24]; % Connection Probability, given by Ko et al 2011
Directional_Probability = [.1 .1; .4 .1]; % Unidirectional 0 degrees = 10%, 90 degrees = 10% ; Bidirectional 0 degrees = 40%, 90 degrees = 10%
Neuron_Connected = []; k =1;
for i = 1:length(neuron.motion.number)
    for j = 1:length(neuron.motion.number)
        if Coding_Probability(neuron.motion.direction(i)) > rand(1) && i~=j
            Neuron_Connected(k,1) = neuron.motion.number(i); % Source Neuron
            Neuron_Connected(k,2) = neuron.motion.number(j); % Connected Neuron
            if Coding_Probability(neuron.motion.direction(i)) == 10/26
                Neuron_Connected(k,3) = (.8 > rand(1)); % 1 = Bidirectional, 0 = Unidirectional in direction of connected pair
            elseif Coding_Probability(neuron.motion.direction(i)) == 4/24
                Neuron_Connected(k,3) = (.5 > rand(1));
            end
            k = k+1;
        end
    end
end

Neuron_Inhibitory_Oscillatory = neuron.inhibitory(randperm(numel(neuron.inhibitory), ceil(length(neuron.inhibitory)*Inhibitory_Oscillatory_Percentage)));
Neuron_Inhibitory_Non_Oscillatory = setdiff(neuron.inhibitory,Neuron_Inhibitory_Oscillatory);
Neuron_Oscillatory_Type = zeros(1,NumNeurons); Neuron_Oscillatory_Type(Neuron_Inhibitory_Oscillatory) = 1; % Stores if neuron is non-oscillatory (0), or oscillatory (1)

%% Motif Location Selector
Neuron_Population_Matrix = zeros(sx, sy); % Initialize location matrix size nxm

ii = 1; iii = 1; MotifEdgeDistance = 50;
for i = 1:NumPads
    
    pad.start.x(i) = 1 + (ii-1).*(sx/NumPadsX); % X start of pad i
    pad.end.x(i) = (ii).*(sx/NumPadsX) - 1; % X end value of pad i
    
    pad.start.y(i) = 1 + (iii-1).*(sy/NumPadsY); % Y start of pad i
    pad.end.y(i) = (iii).*(sy/NumPadsY) - 1; % Y end value of pad i
    
    % modstart & modend filters so that motifs are not placed too close to
    % edge of xy matrix
    pad.modstart.x(i) = pad.start.x(i);
    pad.modend.x(i) = pad.end.x(i);
    pad.modstart.y(i) = pad.start.y(i);
    pad.modend.y(i) = pad.end.y(i);
    
    % X Filtering
    if sum(ii == 1:NumPadsX:NumPads) > 0 % If this is a top pad
        pad.modstart.x(i) = pad.start.x(i) + MotifEdgeDistance;
    elseif sum(ii == NumPadsX:NumPadsX:NumPads) > 0 % if this is a bottom pad
        pad.modend.x(i) = pad.end.x(i) - MotifEdgeDistance;
    end
    
    % Y Filtering
    if sum(ii == 1:NumPadsX) > 0 % If this is a top pad
        pad.modstart.y(i) = pad.start.y(i) + MotifEdgeDistance;
    elseif sum(ii == NumPads-NumPadsX+1:NumPads) > 0 % if this is a bottom pad
        pad.modend.y(i) = pad.end.y(i) - MotifEdgeDistance;
    end
    
    % Next Pad Number
    ii = ii +1;
    if ii > NumPadsX % Once we have reached the last column, reset x and go to the next y row of pads
        ii = 1;
        iii = iii +1;
    end
end

ii = 1; ClosestMotif = 0; MinimumMotifDistance = 100;
for i = 1:NumMotifs
    
    motif.pad(i) = ii;
    motif.orientation(i) = randi([1 4],1,1); % Assigns orientation. 0,90,180,270 degree change
    
    if i == 1 % First pad only, initialize random location
        motif.center.x(i) = randi([pad.modstart.x(ii) pad.modend.x(ii)],1,1);
        motif.center.y(i) = randi([pad.modstart.y(ii) pad.modend.y(ii)],1,1);
    else
        while ClosestMotif < MinimumMotifDistance
            x1 = randi([pad.modstart.x(ii) pad.modend.x(ii)],1,1); % Selects random X,Y Combination to test
            y1 = randi([pad.modstart.y(ii) pad.modend.y(ii)],1,1);
            ClosestMotif = min(sqrt((x1 - motif.center.x).^2 + (y1 - motif.center.y).^2)); % Determines ecleudian distance between test xy and all motifs
            if ClosestMotif > MinimumMotifDistance % If ecleudian distance to every other motif is greater than minimum, these values work.
                motif.center.x(i) = x1;
                motif.center.y(i) = y1;
            end
        end
    end
    
    ClosestMotif = 0; % Reset value for closest motif distance
    
    if ii+1 > NumPads % Evenly distributes motifs across pads. resets pad # if above max
        ii = 1; % Reset pad sequence
    else
        ii = ii + 1; % Next pad in sequence
    end
    
end


%% Neural Population Matrices

MinimumNeuronDistance = neuron.radii*2 + 1;
MotifLength = 50; ClosestNeuron = 0;
for i = 1:NumNeurons
    x1 = motif.center.x(neuron.motif(i));
    y1 = motif.center.y(neuron.motif(i));
    
    if i == 1 % First neuron only, initialize random location
        neuron.x(i) = randi([x1-MotifLength x1+MotifLength],1,1); % Selects random X,Y Combination to test
        neuron.y(i) = randi([y1+neuron.radii+1 y1+MotifLength],1,1);
    else
        while ClosestNeuron < MinimumNeuronDistance
            if neuron.type == 1
                y2 = randi([y1+neuron.radii+1 y1+MotifLength],1,1);
            else
                y2 = randi([y1-MotifLength y1-neuron.radii-1],1,1);
            end
            x2 = randi([x1-MotifLength x1+MotifLength],1,1); % Selects random X,Y Combination to test
            ClosestNeuron = min(sqrt((x2 - neuron.x).^2 + (y2 - neuron.y).^2)); % Determines ecleudian distance between test xy and all motifs
            if ClosestNeuron > MinimumNeuronDistance % If ecleudian distance to every other motif is greater than minimum, these values work.
                neuron.x(i) = x2;
                neuron.y(i) = y2;
            end
        end
    end
    ClosestNeuron = 0; % Reset value for closest motif distance
end

% Create neuron.axon for inhibitory neurons
neuron.number.indices.axon = {};

for i = 1:NumMotifs
    axonmap = zeros(sx,sy); % reset axonmap
    % Axon hub line
    x = find(neuron.motif == i);
    xc = [min(neuron.x(x)) max(neuron.x(x))];
    yc = [motif.center.y(i) motif.center.y(i)];
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind(size(axonmap), rIndex, cIndex); % Linear indices
    axonmap(lindex) = 1;
    
    % Axon connectiotion from neuron to hub line
    for ii = 1:length(x)
        xc = [neuron.x(x(ii)) neuron.x(x(ii))];
        yc = [neuron.y(x(ii)) motif.center.y(i)];
        Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
        rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
        cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
        lindex = sub2ind(size(axonmap), rIndex, cIndex); % Linear indices
        axonmap(lindex) = 1;
    end
    
    motif.axon = find(axonmap == 1); % Creates motif.axon, the inhibitory hub connection
    for ii = 1:length(x)
        if neuron.type(x(ii)) == 1
            neuron.number.indices.axon = motif.axon;
        end
    end
end

%% Stimulation Locations Matrix

electrode.x = []; electrode.y = [];
for i = 1:10 % Number of electrodes = 100, 10x10 over 4mm in center of map. Does not create an electrode at percent center!
    for j = 1:10
        electrode.x(i,j) = (i-1) * sx/10 + sy/20;
        electrode.y(i,j) = (j-1) * sy/10 + sy/20;
    end
end
electrode.x = electrode.x(:); electrode.y = electrode.y(:);

Stim_Distance_Map = zeros(length(electrode.x),sx,sy);
for i = 1:length(electrode.x)
    Stim_Loc =  zeros(sx, sy); % Location of stimulus
    Stim_Loc(electrode.y(i)-Electrode.radii:electrode.y(i)+Electrode.radii,electrode.x(i)-Electrode.radii:electrode.x(i)+Electrode.radii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
    Stim_Distance_Map(i,:,:) = Stim_Distance;
end

%% Neuron And Axon Area Summation

Axon_Excitatory_matrix = zeros(sx,sy);
neuron.number.i0.soma = zeros(NumNeurons,length(electrode.x));
neuron.number.i0.axon = zeros(NumNeurons,length(electrode.x)); % Stores 1/r^2 area component of Axons
neuron.number.direction = randi([1,4],NumNeurons,1); % 1 = x value down (left),2 = Y value down (up), 3 = X value up (right), 4 = Y value up (down)
neuron.number.direction(neuron.inhibitory) = motif.orientation; % Sets stored value for inhibitory axon direction


neuron.number.indices.soma = {};
for i = 1:NumNeurons
    neuron.number(i).indices.soma = {sub2ind([sx sy],NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii)};
end

for i = 1:NumNeurons
    I0_Axon_Neurons1 = zeros(1,length(electrode.x));
    I0_Soma_Neurons1 = zeros(1,length(electrode.x));
    if neuron.type(i) == 1
        x = MotifX(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifX(Neuron_Motif(i))+length(Neuron_Loc)/2-1; % Location of motif does not change with Orientation
        y = MotifY(Neuron_Motif(i))-length(Neuron_Loc)/2:MotifY(Neuron_Motif(i))+length(Neuron_Loc)/2-1;
        lindex = find(Axon_Inhibitory_matrix(y,x)==1);
        for j = 1:length(electrode.x)
            I0_Axon_Neurons1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % sums all inhibitory axon current feedback
            I0_Soma_Neurons1(j) = sum(sum(1./Stim_Distance_Map(j,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma % Summation 1/r^2 area component of soma
        end
        
    elseif neuron.type(i) == 2 % Creating a new axon matrix for Excitatory neurons. Must be made one at a time to calculate current effect on every neuron
        re = [0,0,0,0]; re(neuron.number(i).direction) = 1; % Random element selection to choose orientation of axon
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
        for j = 1:length(electrode.x)
            I0_Axon_Neurons1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
            I0_Soma_Neurons1(j) = sum(sum(1./Stim_Distance_Map(j,NeuronY(i)-NeuronRadii:NeuronY(i)+NeuronRadii, NeuronX(i)-NeuronRadii:NeuronX(i)+NeuronRadii).^2)); % Summation 1/r^2 area component of soma
        end
        
    end
    neuron.number(i).indices.axon = lindex(:);
    neuron.number.i0.axon(i,:) = I0_Axon_Neurons1;
    neuron.number.i0.soma(i,:) = I0_Soma_Neurons1;
end

% Axon Generation for motion-tuned neural connections
I0_Motion_Neuron_Connections = zeros(length(Neuron_Connected),length(electrode.x));
for i = 1:length(Neuron_Connected)
    I0_Motion_Neuron_Connections1 = zeros(1,length(electrode.x));
    xc = [NeuronX(Neuron_Connected(i,1)) NeuronX(Neuron_Connected(i,2))]; % x coordinates
    yc = [NeuronY(Neuron_Connected(i,1)) NeuronY(Neuron_Connected(i,2))]; % y coordinates
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind([sx sy], rIndex,cIndex); % Linear indices
    for j = 1:length(electrode.x)
        I0_Motion_Neuron_Connections1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
    end
    I0_Motion_Neuron_Connections(i,:) = I0_Motion_Neuron_Connections1(j);
end
neuron.number(i).i0.Motion = zeros(NumNeurons,length(electrode.x));
for i = 1:length(I0_Motion_Neuron_Connections)
    neuron.number(Neuron_Connected(i,1)).i0.motion = neuron.number(Neuron_Connected(i,1)).i0.motion + I0_Motion_Neuron_Connections(i,:);
    neuron.number(Neuron_Connected(i,2)).i0.motion = neuron.number(Neuron_Connected(i,2)).i0.motion + I0_Motion_Neuron_Connections(i,:);
end

Axon_Excitatory_matrix(neuron.number(neuron.excitatory).indices.axon) = 1; % Builds excitatory axon matrix


neuron.number.currentmult = zeros(NumNeurons,length(electrode.x));
Axon_Direction_Theta = zeros(NumNeurons,1);
for i = 1:NumNeurons
    if neuron.number(i).direction == 1 %(left = 180 degrees)
        Axon_Direction_Theta(i) = 180;
    elseif neuron.number(i).direction == 2 % (up = 90 degrees)
        Axon_Direction_Theta(i) = 90;
    elseif neuron.number(i).direction == 3 % (right = 0 degrees)
        Axon_Direction_Theta(i) = 0;
    elseif neuron.number(i).direction == 4 % (down, = -90 or 270 degrees)
        Axon_Direction_Theta(i) = -90;
    end
end

for i = 1:NumNeurons
    x2 = NeuronX(i);
    y2 = NeuronY(i);
    
    for ii = 1:length(electrode.x)
        
        x1 = electrode.x(ii);
        y1 = electrode.y(ii);
        theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
        normDeg = mod(Axon_Direction_Theta(i)-theta,360);
        absDiffDeg = min(360-normDeg, normDeg);
        
        if absDiffDeg < theta_threshold % If angle meets critera then:
            neuron.number(i).currentmult(ii) = 1; % Stores directional current multiplier as full
        else % If angle is not within criteria
            neuron.number(i).currentmult(ii) = 0.25; % Stores current summation as depreciated
        end
    end
end

%% Cleaning output

clear ('Stim_Distance_Map','rp','rp_x','rp_y','Stim_Distance','Stim_Loc','PadSize','Neuron_points_Matrix','Ed');
save InitialConditionsFull.mat;
