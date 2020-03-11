clear all;

%% Initial Condition Parameters

% Neuron properties
neuron.radii = 5; % in micrometers
Electroderadii = 10; % in micrometers
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
neuron.number.axon = {};
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
            neuron.number(i).axon = motif.axon;
        end
    end
end



