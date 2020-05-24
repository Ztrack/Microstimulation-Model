clear all; clc;

%% Initial Condition Parameters

% Neuron properties
neuron.radii = 5; % in micrometers
electrode.radii = 10; % in micrometers
Inhibitory_Factor = .01; % 0 to 1 Greater factor will inhibit RS neurons more
theta_threshold = 45; % angle difference threshold - If the neuron axon is out of phase by at least this much to the current-stimulus, the electric field acting on the neuron is weakened to 25%.

% Population Properties
NumNeurons = 1000; % Must be multiples of 5 to satisfy ratio, if using 1:4 ratio. Every 1st neuron is Inhibitory
NumNeuronsMotif = 10; % # of neurons in each motif
NumInhibitoryMotif = 5; % # of inhibitory neurons in motif
NumMotifs = NumNeurons/NumNeuronsMotif; % # Neuron in Motif.
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
NumSteps = 1/dt;

%% Neuron Definitions

count = ones(1,NumNeuronsMotif); count(NumInhibitoryMotif+1:end) = 2;
neuron.type = repmat([count],1,NumMotifs); % Defines what type of neuron this number is. Inhibitory =1, Excitatory = 2
neuron.motif = repmat(1:1:NumMotifs,NumNeuronsMotif,1); neuron.motif = neuron.motif(:)'; % Defines which motif the neuron belongs to
neuron.inhibitory = find(neuron.type == 1); % Array with all inhibitory neurons
neuron.excitatory = find(neuron.type == 2); % Array with all Excitatory neurons
neuron.motion.number = sort(neuron.excitatory(randperm(length(neuron.excitatory),floor(length(neuron.excitatory)*NeuronMotionRatio)))); % Neurons selected for motion
neuron.nonmotion.number = ones(NumNeurons,1); neuron.nonmotion.number(neuron.motion.number) = 0; neuron.nonmotion.number = find(neuron.nonmotion.number == 1); % Non-Motion neurons
neuron.motion.direction = randi([1,2],[1,length(neuron.motion.number)]); % Gives every motion neuron an orientation 0 or 90 degrees


% Coding_Probability = [10/26 4/24]; % Connection Probability, given by Ko et al 2011
% Directional_Probability = [.1 .1; .4 .1]; % Unidirectional 0 degrees = 10%, 90 degrees = 10% ; Bidirectional 0 degrees = 40%, 90 degrees = 10%
% Neuron_Connected = []; k =1;
% for i = 1:length(neuron.motion.number)
%     for j = 1:length(neuron.motion.number)
%         if Coding_Probability(neuron.motion.direction(i)) > rand(1) && i~=j
%             Neuron_Connected(k,1) = neuron.motion.number(i); % Source Neuron
%             Neuron_Connected(k,2) = neuron.motion.number(j); % Connected Neuron
%             if Coding_Probability(neuron.motion.direction(i)) == 10/26
%                 Neuron_Connected(k,3) = (.8 > rand(1)); % 1 = Bidirectional, 0 = Unidirectional in direction of connected pair
%             elseif Coding_Probability(neuron.motion.direction(i)) == 4/24
%                 Neuron_Connected(k,3) = (.5 > rand(1));
%             end
%             k = k+1;
%         end
%     end
% end

neuron.oscillatory = neuron.inhibitory(randperm(numel(neuron.inhibitory), ceil(length(neuron.inhibitory)*Inhibitory_Oscillatory_Percentage)));
neuron.nonoscillatory = setdiff(neuron.inhibitory,neuron.oscillatory);
neuron.oscillatorytype = zeros(1,NumNeurons); neuron.oscillatorytype(neuron.oscillatory) = 1; % Stores if neuron is non-oscillatory (0), or oscillatory (1)

%% Motif Location Selector

% The x,y space is divided into x number of pads. WE must determine where
% these pads lie. For the purposes of keeping things within bounds, a
% modstart and modend value is determined. Motif locations are then
% selected depending on pad designation and distance to other motifs.

% Pad location Definitions
ii = 1; iii = 1; MotifEdgeDistance = 200;
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
    end
    if sum(ii == NumPadsX:NumPadsX:NumPads) > 0 % if this is a bottom pad
        pad.modend.x(i) = pad.end.x(i) - MotifEdgeDistance;
    end
    
    % Y Filtering
    if sum(iii == 1) > 0 % If this is a top pad
        pad.modstart.y(i) = pad.start.y(i) + MotifEdgeDistance;
    end
    if sum(iii == NumPadsY) > 0 % if this is a bottom pad
        pad.modend.y(i) = pad.end.y(i) - MotifEdgeDistance;
    end
    
    % Next Pad Number
    ii = ii +1;
    if ii > NumPadsX % Once we have reached the last column, reset x and go to the next y row of pads
        ii = 1;
        iii = iii +1;
    end
end

% Motif Location Designation
ii = 1; ClosestMotif = 0; MotifLength = 50; MinimumMotifDistance = MotifLength*2;
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
        end
        motif.center.x(i) = x1; % Values are accepted, they are far enough away from all other motifs
        motif.center.y(i) = y1;
    end
    
    ClosestMotif = 0; % Reset value for closest motif distance
    
    if ii+1 > NumPads % Evenly distributes motifs across pads. resets pad # if above max
        ii = 1; % Reset pad sequence
    else
        ii = ii + 1; % Next pad in sequence
    end
    
end


%% Neural Population Matrices

MinimumNeuronDistance = neuron.radii*4 + 1;
ClosestNeuron = 0; 
neuron.x(i) = 1; neuron.y(i) = 1; % Initialize
for i = 1:NumNeurons
    x1 = motif.center.x(neuron.motif(i)); % X Value of motif center
    y1 = motif.center.y(neuron.motif(i)); % Y value of motif center
    
    while ClosestNeuron < MinimumNeuronDistance
        % Selects random X,Y Combination to test
        r = randi([neuron.radii*3 MotifLength],1,1);
        
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
    
    % Axon connection from neuron to hub line
    for ii = 1:length(x)
        xc = [neuron.x(x(ii)) neuron.x(x(ii))]; % Line does not move in x direction
        if (neuron.y(x(ii)) - motif.center.y(i)) > 0 % If the neuron is above the center
            yc = [neuron.y(x(ii))-neuron.radii-1 motif.center.y(i)]; % Create a line down to center
        else
            yc = [neuron.y(x(ii))+neuron.radii+1 motif.center.y(i)]; % Otherwise create line up to center
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

for i = 1:NumMotifs
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
            [y3,x3] = ind2sub([sy sx],neuron.indices(x(ii)).axon); % Converts indices to x,y. returns rows,columns
            x3 = motif.center.x(i) - x3;
            y3 = motif.center.y(i) - y3;
            x2 = motif.center.x(i) + x3*cosd(theta)-y3*sind(theta); % Transform x deg
            y2 = motif.center.y(i) + x3*sind(theta)+y3*cosd(theta); % Transform x deg
            ind = sub2ind([sy sx],y2,x2);
            
            neuron.indices(x(ii)).axon = ind; % Converts x,y back to indices for storage
            
        end
    end
end

% Excitatory Axon Generation
for i = 1:NumNeurons
    if neuron.type(i) == 1
        neuron.direction(i) = motif.orientation(neuron.motif(i)); % Sets stored value for inhibitory axon direction
    elseif neuron.type(i) == 2 % Creating a new axon matrix for Excitatory neurons. Must be made one at a time to calculate current effect on every neuron
        neuron.direction(i) = randi([1,4],1,1); % 1 = x value down (left),2 = Y value down (up), 3 = X value up (right), 4 = Y value up (down)
        re = [0,0,0,0]; re(neuron.direction(i)) = 1; % Random element selection to choose orientation of axon
        rl = randi([50+neuron.radii,480+neuron.radii],1); % Random length of axon + neuron.radii
        xc = [neuron.x(i)+(re(3)*neuron.radii)-(re(1)*neuron.radii) neuron.x(i)+(re(3)*rl)-(re(1)*rl)]; % x coordinates of neuron, axon point
        yc = [neuron.y(i)+(re(4)*neuron.radii)-(re(2)*neuron.radii) neuron.y(i)+(re(4)*rl)-(re(2)*rl)]; % y coordinates of neuron, axon point
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
        neuron.indices(i).axon = lindex(:);
    end
end

% Stores which pad a neuron belongs to

for i = 1:NumNeurons
   neuron.pad(i) = motif.pad(neuron.motif(i));
end

% Neuron Soma Indices

for i = 1:NumNeurons
    singlepop = zeros(sx,sy);
    singlepop(neuron.y(i)-neuron.radii:neuron.y(i)+neuron.radii, neuron.x(i)-neuron.radii:neuron.x(i)+neuron.radii) = 1;
    neuron.indices(i).soma = find(singlepop == 1);
end

clear('singlepop');
%% Stimulation Locations Matrix

% 10x10 electrode array. First we must determine their locations
electrode.x = []; electrode.y = [];
for i = 1:10 % Number of electrodes = 100, 10x10 over 4mm in center of map. Does not create an electrode at percent center!
    for j = 1:10
        electrode.x(i,j) = (i-1) * sx/10 + sy/20;
        electrode.y(i,j) = (j-1) * sy/10 + sy/20;
    end
end
electrode.x = electrode.x(:); electrode.y = electrode.y(:);

% Now we must determine how far each square is to each electrode for later
% intergrating.

Stim_Distance_Map = zeros(length(electrode.x),sx,sy);
for i = 1:length(electrode.x)
    Stim_Loc =  zeros(sx, sy); % Location of stimulus
    Stim_Loc(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1; % Builds electrode stimulus location matrix
    Ed = bwdist(Stim_Loc); % Calculates the Euclidean distance from stim points. Sets origin (stim point) to 0
    Stim_Distance = Stim_Loc + Ed; % Distance of every square to the nearest stim location
    Stim_Distance_Map(i,:,:) = Stim_Distance;
end

%% Neuron And Axon Area Summation
load('lightspread.mat');
neuron.io.soma = zeros(NumNeurons,length(electrode.x)); % Stores 1/r^2 area component of Somas
neuron.io.axon = zeros(NumNeurons,length(electrode.x)); % Stores 1/r^2 area component of Axons
neuron.oo.soma = zeros(NumNeurons,length(electrode.x)); % Stores light stimulus area summation component of Somas
lightspread.calc = [];
for i = 1:NumNeurons

    for j = 1:length(electrode.x)
        
        neuron.io.axon(i,j) = sum(sum(1./Stim_Distance_Map(j,neuron.indices(i).axon).^2)); % Summation 1/r^2 area component of axon
        neuron.io.soma(i,j) = sum(sum(1./Stim_Distance_Map(j,neuron.indices(i).soma).^2)); % Summation 1/r^2 area component of soma % Summation 1/r^2 area component of soma
        
        lightspread.calc = lightspread.averaged.a.*exp(lightspread.averaged.b.*Stim_Distance_Map(j,neuron.indices(i).soma)); % Light Intensity Calculation
        % lightspread.calc(lightspread.calc < 0) = 0; % Should not happen, debugging
        neuron.oo.soma(i,j) = sum(sum(lightspread.calc));
        
    end
    
end

% Axon Generation for motion-tuned neural connections
% I0_Motion_Neuron_Connections = zeros(length(Neuron_Connected),length(electrode.x));
% for i = 1:length(Neuron_Connected)
%     I0_Motion_Neuron_Connections1 = zeros(1,length(electrode.x));
%     xc = [neuron.x(Neuron_Connected(i,1)) neuron.x(Neuron_Connected(i,2))]; % x coordinates
%     yc = [neuron.y(Neuron_Connected(i,1)) neuron.y(Neuron_Connected(i,2))]; % y coordinates
%     Points = max(abs(diff(xc)), abs(diff(yc)))+1; % Number of points in line
%     rIndex = round(linspace(yc(1), yc(2), Points)); % Row indices
%     cIndex = round(linspace(xc(1), xc(2), Points)); % Column indices
%     lindex = sub2ind([sx sy], rIndex,cIndex); % Linear indices
%     for j = 1:length(electrode.x)
%         I0_Motion_Neuron_Connections1(j) = sum(sum(1./(Stim_Distance_Map(j,lindex).^2))); % Stores information on the axon morphology to calculate current backpropogation
%     end
%     I0_Motion_Neuron_Connections(i,:) = I0_Motion_Neuron_Connections1(j);
% end
% neuron.number(i).i0.Motion = zeros(NumNeurons,length(electrode.x));
% for i = 1:length(I0_Motion_Neuron_Connections)
%     neuron.number(Neuron_Connected(i,1)).i0.motion = neuron.number(Neuron_Connected(i,1)).i0.motion + I0_Motion_Neuron_Connections(i,:);
%     neuron.number(Neuron_Connected(i,2)).i0.motion = neuron.number(Neuron_Connected(i,2)).i0.motion + I0_Motion_Neuron_Connections(i,:);
% end

%% Population Maps

population.inhibitory.map = zeros(sx,sy); % Inhibitory neuron population map
population.excitatory.map = zeros(sx,sy); % Excitatory neuron population map
population.axon.map = zeros(sx,sy); % Axon only population map
population.oscillatory.map = zeros(sx,sy); % Population of oscillating neurons
population.motion.map = zeros(sx,sy); % Motion-Tuned Neuron Population
population.all.map = zeros(sx,sy); % population of all elements

for i = 1:NumNeurons
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

neuron.dirmult = zeros(NumNeurons,length(electrode.x));
Axon_Direction_Theta = zeros(NumNeurons,1);
for i = 1:NumNeurons
    if neuron.direction(i) == 1 %(left = 180 degrees)
        Axon_Direction_Theta(i) = 180;
    elseif neuron.direction(i) == 2 % (up = 90 degrees)
        Axon_Direction_Theta(i) = 90;
    elseif neuron.direction(i) == 3 % (right = 0 degrees)
        Axon_Direction_Theta(i) = 0;
    elseif neuron.direction(i) == 4 % (down, = -90 or 270 degrees)
        Axon_Direction_Theta(i) = -90;
    end
end

% Calculate angle from axon hillic to each electrode

for i = 1:NumNeurons
    x2 = neuron.x(i);
    y2 = neuron.y(i);
    
    for ii = 1:length(electrode.x)
        
        x1 = electrode.x(ii);
        y1 = electrode.y(ii);
        theta = atan2d(y2-y1,x2-x1); % Angle from -180 to 180
        normDeg = mod(Axon_Direction_Theta(i)-theta,360);
        absDiffDeg = min(360-normDeg, normDeg);
        
        if absDiffDeg < theta_threshold % If angle meets critera then:
            neuron.dirmult(i,ii) = 1; % Stores directional current multiplier as full
        else % If angle is not within criteria
            neuron.dirmult(i,ii) = 0.25; % Stores current summation as depreciated
        end
    end
end

%% Cleaning output

clear ('Stim_Distance_Map','rp','rp_x','rp_y','Stim_Distance','Stim_Loc','PadSize','Neuron_points_Matrix','Ed');
save InitialConditionsFullB.mat;