clearvars; clc;

%% Initial Condition Parameters

% Neuron properties
params.neuron.radii = 5; % in micrometers
params.electrode.radii = 10; % in micrometers
params.theta_threshold = 90; % angle difference threshold - If the neuron axon is out of phase by at least this much to the current-stimulus, the electric field acting on the neuron is weakened to 25%.
params.axonhillockradii = 1;

% Population Properties
params.numneurons = 1000; % Must be multiples of 5 to satisfy ratio, if using 1:4 ratio. Every 1st neuron is Inhibitory
params.numelectrodes = 100;
params.NumNeuronsMotif = 5; % # of neurons in each motif
params.numinhibitorymotif = 1; % # of inhibitory neurons in motif
params.NumMotifs = params.numneurons/params.NumNeuronsMotif; % # Neuron in Motif.
params.inhibitory_oscillatory_percentage = .33; % Percentage of inhibitory neurons which oscillate
params.NeuronMotionRatio = 0.2; % Ratio of Motion Neurons to non-motion neurons. based on population of excitatory neurons

% Spatial Properties
params.sx = 4800; % Size of work area
params.sy = 4800; % Size of work area
params.numpadsx = 5; % Number of pads in X Direction
params.numpadsy = 3; % Number of pads in Y direction
params.numpads = params.numpadsx*params.numpadsy; % Total number of pads in params.sx*params.sy area

% Temporal Properties
T = 1; % Duration of experiment in seconds
dt = T/1000; % Time step length
numsteps = 1/dt;

%% Neuron Definitions

count = ones(1,params.NumNeuronsMotif); count(params.numinhibitorymotif+1:end) = 2;
neuron.type = repmat([count],1,params.NumMotifs); % Defines what type of neuron this number is. Inhibitory =1, Excitatory = 2
neuron.motif = repmat(1:1:params.NumMotifs,params.NumNeuronsMotif,1); neuron.motif = neuron.motif(:)'; % Defines which motif the neuron belongs to
neuron.inhibitory = find(neuron.type == 1); % Array with all inhibitory neurons
neuron.excitatory = find(neuron.type == 2); % Array with all Excitatory neurons
neuron.motion.number = sort(neuron.excitatory(randperm(length(neuron.excitatory),floor(length(neuron.excitatory)*params.NeuronMotionRatio)))); % Neurons selected for motion
neuron.nonmotion.number = ones(params.numneurons,1); neuron.nonmotion.number(neuron.motion.number) = 0; neuron.nonmotion.number = find(neuron.nonmotion.number == 1); % Non-Motion neurons
neuron.motion.tuning = randi([1 8],size(neuron.motion.number));

% Neuron Oscillation Properties
neuron.oscillatory = neuron.inhibitory(randperm(numel(neuron.inhibitory), ceil(length(neuron.inhibitory)*params.inhibitory_oscillatory_percentage)));
neuron.nonoscillatory = setdiff(neuron.inhibitory,neuron.oscillatory);
neuron.oscillatorytype = zeros(1,params.numneurons); neuron.oscillatorytype(neuron.oscillatory) = 1; % Stores if neuron is non-oscillatory (0), or oscillatory (1)

% Neuron Adaptation
neuron.adapt.type = ones(1,params.numneurons); % 1 = Fast Adapting Neurons
randorder = randperm(length(neuron.adapt.type)); % Random Permutation
neuron.adapt.type(randorder(1:floor(params.numneurons.*.06))) = 2; % 2 = Slow Adapting Neurons
randorder = setdiff(randperm(length(neuron.adapt.type)),find(neuron.adapt.type == 2)); % Random Permutation
neuron.adapt.type(randorder(1:floor(params.numneurons.*.40))) = 3; % Mixed
load('rateFunctions.mat');
neuron.adapt.ratefunction(1,:) = RA_function;
neuron.adapt.ratefunction(2,:) = SA_function;
neuron.adapt.ratefunction(3,:) = MIXED_function;

%% Motif Location Selector

% The x,y space is divided into x number of pads. WE must determine where
% these pads lie. For the purposes of keeping things within bounds, a
% modstart and modend value is determined. Motif locations are then
% selected depending on pad designation and distance to other motifs.

% Pad location Definitions
pad.PadsXLength = params.sx/params.numpadsx; % Length of a X pad
pad.PadsYLength = params.sy/params.numpadsy; % Length of a Y pad
PadSize = zeros(pad.PadsYLength, pad.PadsXLength);
pad.center.x = repmat(1:params.numpadsx,1,params.numpadsy)*pad.PadsXLength;
pad.center.y = sort(repmat(1:params.numpadsy,1,params.numpadsx))*pad.PadsYLength;
ii = 1; iii = 1; MotifEdgeDistance = 200;
for i = 1:params.numpads
    
    pad.start.x(i) = 1 + (ii-1).*(params.sx/params.numpadsx); % X start of pad i
    pad.end.x(i) = (ii).*(params.sx/params.numpadsx) - 1; % X end value of pad i
    
    pad.start.y(i) = 1 + (iii-1).*(params.sy/params.numpadsy); % Y start of pad i
    pad.end.y(i) = (iii).*(params.sy/params.numpadsy) - 1; % Y end value of pad i
    
    % modstart & modend filters so that motifs are not placed too close to
    % edge of xy matrix
    pad.modstart.x(i) = pad.start.x(i);
    pad.modend.x(i) = pad.end.x(i);
    pad.modstart.y(i) = pad.start.y(i);
    pad.modend.y(i) = pad.end.y(i);
    
    % X Filtering
    if sum(ii == 1:params.numpadsx:params.numpads) > 0 % If this is a top pad
        pad.modstart.x(i) = pad.start.x(i) + MotifEdgeDistance;
    end
    if sum(ii == params.numpadsx:params.numpadsx:params.numpads) > 0 % if this is a bottom pad
        pad.modend.x(i) = pad.end.x(i) - MotifEdgeDistance;
    end
    
    % Y Filtering
    if sum(iii == 1) > 0 % If this is a top pad
        pad.modstart.y(i) = pad.start.y(i) + MotifEdgeDistance;
    end
    if sum(iii == params.numpadsy) > 0 % if this is a bottom pad
        pad.modend.y(i) = pad.end.y(i) - MotifEdgeDistance;
    end
    
    % Next Pad Number
    ii = ii +1;
    if ii > params.numpadsx % Once we have reached the last column, reset x and go to the next y row of pads
        ii = 1;
        iii = iii +1;
    end
end

% Motif Location Designation
ii = 1; ClosestMotif = 0; MotifLength = 50; MinimumMotifDistance = MotifLength*2;
for i = 1:params.NumMotifs
    
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
        end
        motif.center.x(i) = x1; % Values are accepted, they are far enough away from all other motifs
        motif.center.y(i) = y1;
    end
    
    ClosestMotif = 0; % Reset value for closest motif distance
    
    if ii+1 > params.numpads % Evenly distributes motifs across pads. resets pad # if above max
        ii = 1; % Reset pad sequence
    else
        ii = ii + 1; % Next pad in sequence
    end
    
end


%% Neural Population Matrices

MinimumNeuronDistance = params.neuron.radii*4 + 1;
ClosestNeuron = 0; 
neuron.x(i) = 1; neuron.y(i) = 1; % Initialize
for i = 1:params.numneurons
    x1 = motif.center.x(neuron.motif(i)); % X Value of motif center
    y1 = motif.center.y(neuron.motif(i)); % Y value of motif center
    
    while ClosestNeuron < MinimumNeuronDistance
        % Selects random X,Y Combination to test
        r = randi([params.neuron.radii*3 MotifLength],1,1);
        
        if neuron.type(i) == 1 % If this is an inhibitory neuron
            y2 = y1+r; % it will lie somewhere above the center 
        else
            y2 = y1-r; % it will lie somewhere below the center 
        end
        
        x2 = randi([x1-MotifLength x1+MotifLength],1,1); % Selects random X,Y Combination to test
        ClosestNeuron = min(sqrt((x2 - neuron.x).^2 + (y2 - neuron.y).^2)); % Determines ecleudian distance between test xy and all motifs
    end
    neuron.x(i) = x2; % Neuron values are accepted, they are far enough away from all other neurons
    neuron.y(i) = y2;
    ClosestNeuron = 0; % Reset value for closest motif distance
end

% Create neuron axon indices for inhibitory neurons. First we draw a 'hub'
% line seperating inhibitory and excitatory neurons. This is done for
% visuals mainly. Then we connect the neurons to the hub lines. This 'tree'
% axon network counts as the inhibitory neuron's axon.

for i = 1:params.NumMotifs
    axonmap = zeros(params.sx,params.sy); % reset axonmap
    
    % Axon hub line
    x = find(neuron.motif == i);
    xc = [min(neuron.x(x)) max(neuron.x(x))];
    yc = [motif.center.y(i) motif.center.y(i)];
    Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
    rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
    cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
    lindex = sub2ind(size(axonmap), rIndex, cIndex); % Linear indices
    axonmap(lindex) = 1;
    
    % Axon connection from neuron to hub line
    for ii = 1:length(x)
        xc = [neuron.x(x(ii)) neuron.x(x(ii))]; % Line does not move in x direction
        if (neuron.y(x(ii)) - motif.center.y(i)) > 0 % If the neuron is above the center
            yc = [neuron.y(x(ii))-params.neuron.radii-1 motif.center.y(i)]; % Create a line down to center
        else
            yc = [neuron.y(x(ii))+params.neuron.radii+1 motif.center.y(i)]; % Otherwise create line up to center
        end
        Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
        rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
        cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
        lindex = sub2ind(size(axonmap), rIndex, cIndex); % Linear indices
        axonmap(lindex) = 1;
    end
    
    % Store axon indices
    for ii = 1:length(x)
        if neuron.type(x(ii)) == 1
            neuron.indices(x(ii)).axon = find(axonmap == 1);
        end
    end
end

% Rotate motifs some radom direction using rotation matrix

for i = 1:params.NumMotifs
    x = find(neuron.motif == i); % Neurons within this motif
    theta = motif.orientation(i)*90; % Degree we are rotating by
    for ii = 1:length(x)
        % For neuron points, we need to calculate the x,y with respect to
        % the centerpoint
        x1 = motif.center.x(i)-neuron.x(x(ii));
        y1 = motif.center.y(i)-neuron.y(x(ii));
        x2 = motif.center.x(i) + x1*cosd(theta)-y1*sind(theta); % Transform x deg
        y2 = motif.center.y(i) + x1*sind(theta)+y1*cosd(theta); % Transform x deg
        
        neuron.x(x(ii)) = x2;
        neuron.y(x(ii)) = y2;
        
        if neuron.type(x(ii)) == 1 % Now we rotate the axon indices
            [y3,x3] = ind2sub([params.sy params.sx],neuron.indices(x(ii)).axon); % Converts indices to x,y. returns rows,columns
            x3 = motif.center.x(i) - x3;
            y3 = motif.center.y(i) - y3;
            x2 = motif.center.x(i) + x3*cosd(theta)-y3*sind(theta); % Transform x deg
            y2 = motif.center.y(i) + x3*sind(theta)+y3*cosd(theta); % Transform x deg
            ind = sub2ind([params.sy params.sx],y2,x2);
            
            neuron.indices(x(ii)).axon = ind; % Converts x,y back to indices for storage
            
        end
    end
end

% Excitatory Axon Generation
for i = 1:params.numneurons

    if neuron.type(i) == 1
        neuron.direction(i) = motif.orientation(neuron.motif(i)); % Sets stored value for inhibitory axon direction
        re = [0,0,0,0]; re(neuron.direction(i)) = 1;
        
    elseif neuron.type(i) == 2 % Creating a new axon matrix for Excitatory neurons. Must be made one at a time to calculate current effect on every neuron
        neuron.direction(i) = randi([1,4],1,1); % 1 = x value down (left),2 = Y value down (down), 3 = X value up (right), 4 = Y value up (up)
        re = [0,0,0,0]; re(neuron.direction(i)) = 1; % Random element selection to choose orientation of axon
        rl = randi([50+params.neuron.radii,480+params.neuron.radii],1); % Random length of axon + params.neuron.radii
        
        % Determine Axon start/end locations
        xc = [neuron.x(i)+(re(3)*params.neuron.radii)-(re(1)*params.neuron.radii) neuron.x(i)+(re(3)*rl)-(re(1)*rl)]; % x coordinates of neuron, axon hillock
        yc = [neuron.y(i)+(re(4)*params.neuron.radii)-(re(2)*params.neuron.radii) neuron.y(i)+(re(4)*rl)-(re(2)*rl)]; % y coordinates of neuron, axon hillock
        if xc(2) <= 0 % Simple filter Makes sure the axon end does not go out of bounds in x
            xc(2) = 1;
        elseif xc(2) >= params.sx
            xc(2) = params.sx;
        end
        if yc(2) <= 0 %Simple filter Makes sure the axon end does not go out of bounds in y
            yc(2) = 1;
        elseif yc(2) >= params.sy
            yc(2) = params.sy;
        end
        Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
        rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
        cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
        lindex = sub2ind([params.sx params.sy], rIndex,cIndex); % Linear indices
        neuron.indices(i).axon = lindex(:);
    end
end

% Stores which pad a neuron belongs to

for i = 1:params.numneurons
   neuron.pad(i) = motif.pad(neuron.motif(i));
end

% Neuron Soma Indices

for i = 1:params.numneurons
    singlepop = zeros(params.sx,params.sy);
    singlepop(neuron.y(i)-params.neuron.radii:neuron.y(i)+params.neuron.radii, neuron.x(i)-params.neuron.radii:neuron.x(i)+params.neuron.radii) = 1;
    neuron.indices(i).soma = find(singlepop == 1);
end

% Determine axon direction angle
Axon_Direction_Theta = zeros(params.numneurons,1);
for i = 1:params.numneurons
    if neuron.direction(i) == 1 %(left = 180 degrees)
        Axon_Direction_Theta(i) = 180;
    elseif neuron.direction(i) == 2 % (up = 90 degrees)
        Axon_Direction_Theta(i) = 270;
    elseif neuron.direction(i) == 3 % (right = 0 or 360 degrees)
        Axon_Direction_Theta(i) = 0;
    elseif neuron.direction(i) == 4 % (down, = -90 or 270 degrees)
        Axon_Direction_Theta(i) = 90;
    end
end

% Neuron Membrane Indices
for i = 1:params.numneurons
    
    % Determine membrane locations
    theta = linspace(0, 360, 4*pi*params.neuron.radii);
    x = neuron.x(i) + params.neuron.radii * cosd(theta);
    y = neuron.y(i) + params.neuron.radii * sind(theta);
    xy = round([x', y']);
    xy = unique(xy, 'rows');
    ind = zeros(length(xy),1); % initialize
    for ii = 1:length(xy)
        ind(ii) = sub2ind([params.sy params.sx],xy(ii,1),xy(ii,2));
    end
    neuron.indices(i).membrane = ind; % Store membrane indices
    
    % Determine Axon Hillock
    neuron.direction(i); % 1 = 180deg, 2=90deg, 3=0deg,4=270deg
    xy2 = [neuron.x(i) + params.neuron.radii * cosd(Axon_Direction_Theta(i)) neuron.y(i) + params.neuron.radii * sind(Axon_Direction_Theta(i))];
    p1 = find(xy(:,1) == xy2(1) & xy(:,2) == xy2(2)); % center point of axon hillock on membrane
    p2 = zeros(1+params.axonhillockradii*2,2); % Initialize axon hillock point storage
    p2(1,:) = xy(p1,:);
    % find distance to p1
    %select x closest to p1 as radius
    dist = [];
    for ii = 1:length(xy)
        dist(ii) = sqrt((p2(1,1)-xy(ii,1))^2 + (p2(1,2)-xy(ii,2))^2); % Distance to center of axon hillock
    end
    dist(dist==0) = nan;
    k = 2;
    for ii = 1:params.axonhillockradii
        p1 = find(dist == min(dist));
        p2(k,:) = xy(p1(1),:); k = k+1;
        p2(k,:) = xy(p1(2),:); k = k+1;
    end
    
    ind2 = zeros(length(p2),1); % Initialize
    for ii = 1:length(p2)
        ind2(ii) = sub2ind([params.sy params.sx],xy(ii,1),xy(ii,2));
    end
    neuron.indices(i).axonhillock = ind2;
    
    % Clean membrane indices of axon hillock indices
    neuron.indices(i).membrane = setdiff(neuron.indices(i).membrane, neuron.indices(i).axonhillock);
end

clear('singlepop');
%% Stimulation Locations Matrix

% 10x10 electrode array. First we must determine their locations
electrode.x = []; electrode.y = [];
% for i = 1:10 % Number of electrodes = 100, 10x10 over 4mm in center of map. Does not create an electrode at percent center!
%     for j = 1:10
%         electrode.x(i,j) = (i-1) * params.sx/10 + params.sy/20;
%         electrode.y(i,j) = (j-1) * params.sy/10 + params.sy/20;
%     end
% end
% electrode.x = electrode.x(:); electrode.y = electrode.y(:);

electrode.x = round((params.sx/((sqrt(params.numelectrodes))+1)).*(1:(sqrt(params.numelectrodes)))); electrode.x = repmat(electrode.x,1,sqrt(params.numelectrodes));
electrode.y = round((params.sy/((sqrt(params.numelectrodes))+1)).*(1:(sqrt(params.numelectrodes)))); electrode.y = repmat(electrode.y,sqrt(params.numelectrodes),1); electrode.y = electrode.y(:);

% Now we must determine how far each square is to each electrode for later
% intergrating.
% I/R^2, in units nA/mm^2
% https://journals.physiology.org/doi/full/10.1152/jn.00126.2006

Stim_Distance_Map = zeros(length(electrode.x),params.sx,params.sy);
for i = 1:length(electrode.x)
    Stim_Loc =  zeros(params.sx, params.sy); % Location of stimulus
    Stim_Loc(electrode.y(i)-params.electrode.radii:electrode.y(i)+params.electrode.radii,electrode.x(i)-params.electrode.radii:electrode.x(i)+params.electrode.radii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location in mm
    Stim_Distance_Map(i,:,:) = Stim_Distance./1000; % Stores distance from electrode in mm
end

%% Neuron And Axon Area Summation
load('lightspread.mat');
neuron.io.soma = zeros(params.numneurons,length(electrode.x)); % Stores 1/r^2 area component of Somas
neuron.io.axon = zeros(params.numneurons,length(electrode.x)); % Stores 1/r^2 area component of Axons
neuron.oo.soma = zeros(params.numneurons,length(electrode.x)); % Stores light stimulus area summation component of Somas

lightspread.calc = [];
for i = 1:params.numneurons

    for j = 1:length(electrode.x)
        
        x1a = Stim_Distance_Map(j,neuron.indices(i).axonhillock); % neuron axon hillock Distance Values in mm
        x1b = Stim_Distance_Map(j,neuron.indices(i).membrane); % neuron membrane Distance Values in mm
        x2 = Stim_Distance_Map(j,neuron.indices(i).axon); % Axon Distance Values in mm

        % Sanity check : r = 1.01:0.01:5; plot(r,1./r.^2); % Electrical
        % field as a function of distance
        neuron.io.soma(i,j) = sum(sum(1./(x1a+1).^2))/length(x1a) + (sum(sum(1./(x1b+1).^2))/length(x1b))/length(x1b); % Mean 1/(r+1)^2 area component of cell. Units nA/mm^2
        neuron.io.axon(i,j) = sum(sum(1./(x2+1).^2))/length(x2); % Mean 1/(r+1)^2 area component of axon. Units nA/mm^2
        
        lightspread.calc = lightspread.averaged.a.*(([x1a x1b]*1000).^lightspread.averaged.b); % Light Intensity Calculation in um
        lightspread.calc(lightspread.calc < 0) = 0; % Should not happen, debugging
        neuron.oo.soma(i,j) = sum(sum(lightspread.calc));
        
    end
    
end

%% Neuron Connectivity

% Distance-based connectivity
%The first step is to calculate connection probablility based on pad
%proximity, not actual distance folliwng P = C.*exp(-r); where C is
%constant value and r represents the pad distance
C = 4/24; % Connectivity constant. Within pads at r=0, the probability of E-E connection is C
% Connections are uni-synaptic currently

for i = 1:params.numpads
    for ii = 1:params.numpads
        pad.distance(i,ii) = sqrt(abs(pad.center.x(i)-pad.center.x(ii))^2 + abs(pad.center.y(i)-pad.center.y(ii))^2); % Finds the distance between pads for later use
    end
end
pad.rawdistance = unique(pad.distance);
pad.referencedistance = pad.rawdistance/pad.rawdistance(2); % Divides by the length of a pad 1 distance away for reference

probability_matrix = zeros(params.numneurons,params.numneurons);
for i = 1:params.numneurons
    for ii = 1:params.numneurons
        if i == ii | neuron.type(i) == 1 | neuron.type(ii) == 1
            % if itself or inhibitory do nothing, no connectivity
        else
            paddistance = sqrt(abs(pad.center.x(neuron.pad(i))-pad.center.x(neuron.pad(ii)))^2 + abs(pad.center.y(neuron.pad(i))-pad.center.y(neuron.pad(ii)))^2)/pad.rawdistance(2);
            probability_matrix(i,ii) = C.*exp(-paddistance);
        end
    end
end

%neuron.connectivity defined as (Neuron acted upon,Neuron applying weight)
neuron.connectivity = probability_matrix > rand(size(probability_matrix)); % Stores connectivity data

% Motion-tuned Inhibitory neuron connectivity
%Next we calculate connection probability based on same-feature tuning
neuron.motion.codingprobability = [10/26 11/44 4/24] ;%Connection Probability, given by Ko et al 2011 https://www.nature.com/articles/nature09880
neuron.motion.codingprobability = [neuron.motion.codingprobability zeros(1,length(unique(neuron.motion.tuning))-length(neuron.motion.codingprobability))];

for i = 1:length(neuron.motion.number)
    for ii = 1:length(neuron.motion.number)
        motion_activation = neuron.motion.codingprobability(abs(neuron.motion.tuning(i) - neuron.motion.tuning(ii))+1) > rand(1,1);
        if motion_activation==1 & i~=ii % If neuron activated 
            neuron.connectivity(neuron.motion.number(i),neuron.motion.number(ii)) = 1; % Update connectivity matrix for this combination
        end
    end
end
% NOTE: ADD BICONNECTIVITY STATISTICS HERE

% Inhibitory neuron connectivity % NOTE: CHANGE IF ADDING MORE THAN 1 INHIB
for i = 1:length(neuron.inhibitory)
    neuron.connectivity(neuron.inhibitory(i)+1:neuron.inhibitory(i)+params.NumNeuronsMotif-1,neuron.inhibitory(i)) = 1; % E-I Acting ONTO excitatory neurons
    neuron.connectivity(neuron.inhibitory(i),neuron.inhibitory(i)+1:neuron.inhibitory(i)+params.NumNeuronsMotif-1) = 1; % I-E Acting ONTO this inhibitory neuron
end

%Synaptic Weight values
% The weights (Wji) of these connections were dependent on the type (excitatory or inhibitory) of the presynaptic (i) and postsynaptic (j) neuron
neuron.weight.WIE = 0; % Inhibitory-Excitatory synaptic weight
neuron.weight.WEE = 1/((72.36-41.40)/0.65); % Excitatory-Excitatory synaptic weight
neuron.weight.WEI = -1/((72.36-41.40)/3.48); % Excitatory - Inhibitory synaptic weight 
neuron.weight.WII = 0; % Inhibitory-Inhibitory synaptic weight;

% Now we connect the synaptic weight to each neuron connectivity
neuron.weight.matrix = zeros(params.numneurons,params.numneurons);
for i = 1:params.numneurons
    for ii = 1:params.numneurons
        if neuron.connectivity(i,ii) == 1 % If there is a connection, continue
            if neuron.type(i) == 1 & neuron.type(ii) == 1% If postsynaptic neuron is inhibitory and pre-synaptic is inhibitory, apply WII
                neuron.weight.matrix(i,ii) = neuron.weight.WII;
            elseif neuron.type(i) == 1 & neuron.type(ii) == 2 % Inhibitory - Excitatory
                neuron.weight.matrix(i,ii) = neuron.weight.WIE;
            elseif neuron.type(i) == 2 & neuron.type(ii) == 1 % Excitatory - Inhibitory
                neuron.weight.matrix(i,ii) = neuron.weight.WEI;
            elseif neuron.type(i) == 2 & neuron.type(ii) == 2 % Excitatory - Excitatory
                neuron.weight.matrix(i,ii) = neuron.weight.WEE;
            else
                % else neuron weight = 0
            end
        end
    end
end


%% Population Maps

population.inhibitory.map = zeros(params.sx,params.sy); % Inhibitory neuron population map
population.excitatory.map = zeros(params.sx,params.sy); % Excitatory neuron population map
population.axon.map = zeros(params.sx,params.sy); % Axon only population map
population.oscillatory.map = zeros(params.sx,params.sy); % Population of oscillating neurons
population.motion.map = zeros(params.sx,params.sy); % Motion-Tuned Neuron Population
population.all.map = zeros(params.sx,params.sy); % population of all elements

for i = 1:params.numneurons
    if neuron.type(i) == 1 % Inhibitory Neuron
        if neuron.oscillatorytype(i) == 1 % If true, this neuron oscillates
            population.oscillatory.map(neuron.indices(i).soma) = 1;
        end
        population.inhibitory.map(neuron.indices(i).soma) = 1;
    else % Excitatory Neuron
        if length(intersect(i,neuron.motion.number)) == 1 % If true, this neuron is a motion neuron
            population.motion.map(neuron.indices(i).soma) = 1;
        end
        population.excitatory.map(neuron.indices(i).soma) = 1;
    end
    population.axon.map(neuron.indices(i).axon) = 1; % Creates axon indices population map
end
population.all.map = population.excitatory.map + population.inhibitory.map + population.axon.map;
population.all.map(population.all.map > 1) = 1; % Filters down any values to 1

% Storing Indices for easy plotting
population.inhibitory.indices = find(population.inhibitory.map > 0);
population.excitatory.indices = find(population.excitatory.map > 0);
population.axon.indices = find(population.axon.map > 0);
population.oscillatory.indices = find(population.oscillatory.map > 0);
population.motion.indices = find(population.motion.map > 0);

% Plotting test:
%figure; imagesc(population.all.map); colorbar;

%% Directional Current Multiplier

% The direction of the neuron axon hillic is important. When the axon
% hillic is in the direction of the current stimulus, it is easier for it
% to become excited.

neuron.dirmult = zeros(params.numneurons,length(electrode.x));

% Calculate angle from axon hillic to each electrode

for i = 1:params.numneurons
    x2 = neuron.x(i);
    y2 = neuron.y(i);
    
    for ii = 1:length(electrode.x)
        
        x1 = electrode.x(ii);
        y1 = electrode.y(ii);
        theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
        normDeg = mod(Axon_Direction_Theta(i)-theta,360);
        absDiffDeg = min(360-normDeg, normDeg);
        
        if absDiffDeg < params.theta_threshold % If angle meets critera then:
            neuron.dirmult(i,ii) = 1; % Stores directional current multiplier as full
        else % If angle is not within criteria
            neuron.dirmult(i,ii) = 0.25; % Stores current summation as depreciated
        end
    end
end

%% Cleaning output

clear ('Stim_Distance_Map','rp','rp_x','rp_y','Stim_Distance','Stim_Loc','PadSize','Neuron_points_Matrix','Ed','probability_matrix');
save InitialConditionsFull.mat;
