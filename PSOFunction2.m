function out = PSOFunction(problem, params)

    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function

    nVar = problem.nVar;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables

    %% Parameters of PSO

    MaxIt = params.MaxIt;   % Maximum Number of Iterations
    AdaptiveItMax = params.AdaptiveItMax; % Max number of adaptive iterations
    AdaptiveItThreshold = params.AdaptiveItThreshold; % if the absolute value of it - (it-1) is less than this, we say there is little enough change
    
    nPop = params.nPop;     % Population Size (Swarm Size)

    w = params.w;           % Intertia Coefficient
    wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
    c1 = params.c1;         % Personal Acceleration Coefficient
    c2 = params.c2;         % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    
    
    for ii = 1:length(VarMin)
        MaxVelocity(ii) = 0.2*(VarMax(ii)-VarMin(ii));
        MinVelocity(ii) = -MaxVelocity(ii);
    end
    
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Position = zeros(1,nVar);
    GlobalBest.Cost = inf;

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        %particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
        for ii = 1:length(VarMin)
            particle(i).Position(ii) = unifrnd(VarMin(ii), VarMax(ii),1,1);
        end
        
        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position); % Puts particle positions (electrode AU) through cost function

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position; 
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, 1);


    %% Main Loop of PSO
    
    AdaptiveIt = 0; % If the last ~5 or however many solutions have been the same, we reached a local minima
    for it=1:MaxIt

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            for ii = 1:length(VarMin) % Ensures that positions are within thresholding limits
                particle(i).Velocity(ii) = max(particle(i).Velocity(ii), MinVelocity(ii));
                particle(i).Velocity(ii) = min(particle(i).Velocity(ii), MaxVelocity(ii));
            end
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            for ii = 1:length(VarMin) % Ensures that positions are within thresholding limits
                particle(i).Position(ii) = max(particle(i).Position(ii), VarMin(ii));
                particle(i).Position(ii) = min(particle(i).Position(ii), VarMax(ii));
            end
            
            % Evaluation
            particle(i).Cost = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(it) = GlobalBest.Cost;

        % Display Iteration Information
        if ShowIterInfo
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
        end

        % Damping Inertia Coefficient
        w = w * wdamp;
        
        if it > 1
        if abs(BestCosts(it)-BestCosts(it-1)) < AdaptiveItThreshold
            AdaptiveIt = AdaptiveIt+1;
        else
            AdaptiveIt = 0; % Reset adaptive it
        end
        end
        
        if AdaptiveIt == AdaptiveItMax
            break
        end
    end
    
    out.pop = particle;
    out.BestSol = GlobalBest;
    out.BestCosts = BestCosts;
    out.NumIt = it;
end