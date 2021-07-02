clc 
clearvars
load('InitialConditionsFull.mat')
set(groot,'defaultLineLineWidth',4.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)
%%
lambdatype = 1;
neuron.lambda = zeros(params.numneurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons
NumTrials = 50;
neuron.IC = 0;
CI = 95;
parfor i = 1:NumTrials
    lambdamod(i,:) = IntegrateAndFire(neuron,params,zeros(params.numneurons,1)); % FR of neurons before any stimulation, includes synaptic connections
end
neuron.lambdamod = mean(lambdamod,1)';
%% Single Run
tic

%Setup
rdpso.n = 10; % Size of the population, # of tries or solution sets of 100 current values + 100 opto values
rdpso.nruns=1; % rdpso.nruns is the required number of runs.
rdpso.niter = 500; % Maximum number of iterations.
rdpso.currentub = 10;
rdpso.optoub = 0.01;
lambdatype = 1;

if any(lambdatype==1 | lambdatype==2)
    rdpso.ndim = 100; % Dimension of the problem, # electrodes if lambdatype is 1 or 2
    rdpso.lb = zeros(1,100); %lower limits.
    rdpso.ub = zeros(1,100)+100; %upper limits.
else
    rdpso.ndim = 200; % Dimension of the problem, # electrodes if lambdatype is 1 or 2
    rdpso.lb = zeros(1,200); %lower limits.
    rdpso.ub = zeros(1,200)+100; %upper limits.
end

[sol] = rdpsofun(rdpso,neuron,params,lambdatype);
disp(sol)

figure;
y = mean(sol.costplot);
yerr = 1.96.*std(sol.costplot)/sqrt(length(y));
errorbar(y,yerr);
title('Optimization Function');
xlabel('Iterations');
ylabel('Performance');

toc
%% Multi Run

%Setup
rdpso.n = 100; % Size of the population, # of tries or solution sets of 100 current values + 100 opto values
rdpso.nruns=3; % rdpso.nruns is the required number of runs.
rdpso.niter = 5000; % Maximum number of iterations.
rdpso.currentub = 10;
rdpso.optoub = 0.01;

for i = 1:4
    
    % Apply Setup
    lambdatype = i;
    if i == 1
        rdpso.ndim = 100; % Dimension of the problem, # electrodes x1
        rdpso.lb = zeros(1,100); %lower limits.
        rdpso.ub = zeros(1,100)+100; %upper limits.
        
    elseif i == 2
        rdpso.ndim = 100; % Dimension of the problem, # electrodes x1
        rdpso.lb = zeros(1,100); %lower limits.
        rdpso.ub = zeros(1,100)+100; %upper limits.
        
    elseif any(i==3 | i==4)
        rdpso.ndim = 200; % Dimension of the problem, # electrodes x2
        rdpso.lb = zeros(1,200); %lower limits.
        rdpso.ub = zeros(1,200)+100; %upper limits.
%         rdpso.ub(1:size(rdpso.ub,2)/2) = rdpso.currentub; % Max value for electrical microstimulation
%         rdpso.ub(length(rdpso.ub)/2+1:length(rdpso.ub)) = rdpso.optoub;
    end
        
    % Run algorithm
    sol = rdpsofun(rdpso,neuron,params,lambdatype);
    rdpsosol(i) = sol;
    
end

figure;
for i = 1:4
    y = mean(rdpsosol(i).costplot);
    yerr = 1.96.*std(rdpsosol(i).costplot)/sqrt(length(y));
    errorbar(y,yerr); hold on;
end
title('Optimization Performance');
xlabel('Iteration Number'); ylabel('Nonmotion/Motion Ratio');
legend('MS','Opto','MS+Opto(Excitatory)','MS+Opto(Inhibitory)');
%% Functions

function [rdpsosol] = rdpsofun(rdpso,neuron,params,lambdatype)
% RDPSO code
% This code is written by Wael T. El Sayed, last update was on 12
% March 2018
% %-------------------------------------------------------------------------------------
% Please, If you will use this code, then cite the following papers:
%--------------
% 1- Paper no.1 
%--------------
% Wael T. Elsayed, Yasser G. Hegazy, Mohamed S. El-bages, and Fahmy M. Bendary 
% �Improved Random Drift Particle Swarm Optimization With Self-Adaptive 
% Mechanism for Solving the Power Economic Dispatch Problem,� 
% IEEE Transactions on Industrial Informatics, vol. 13, no. 3, p. 1017 - 1026, 2017.
%-------------
% 2- Paper no.2
%-------------
% J. Sun, V. Palade, X.-J. Wu, W. Fang, and Z. Wang, �Solving the power
% economic dispatch problem with generator constraints by random drift
% particle swarm optimization,� IEEE Transactions on Industrial Informatics,
% vol. 10, no. 1, pp. 222�232, 2014.

% initialize the limits here.

% initialize the popoulation size.

%initialize the algorithm control parameters.
% The range of thermal coefficient parameter see the paper no. 2 above page 226
% second column
alphamax=0.9;
alphamin=0.3;
% Acceleration coefficients
c1=2; 
c2=2; 
% Initialization to the results matrix.
% results=[global cost per run,x1,x2,...,xn, time per run]
results=zeros(rdpso.nruns,rdpso.ndim+2);
rdpsosol.costplot=zeros(rdpso.nruns,rdpso.niter);
jn=0; % Initialization to the loop variable for the runs.

while jn<rdpso.nruns %loop for the number of runs
    jn=jn+1;
    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the population within the limits for the current run.
    pop=zeros(rdpso.n,rdpso.ndim);
    for i=1:rdpso.ndim
    pop(:,i)=rdpso.lb(i)+(rdpso.ub(i)-rdpso.lb(i))*rand(rdpso.n,1);
    end
    
    % Evaluate the total population at the first iteration.
    parfor i = 1:rdpso.n
        [cost(1,i)]=costfunction(pop(i,:),neuron,params,lambdatype,rdpso); % Results in number of pop cost values
    end
    
    %Initialize the local best equal to the initial population.
    localsolution = pop; % location of the local best.
    localcost = cost; % cost of local best.
    
    % Initialize the global solution.
    [globalcost,indx] = min(cost); % Store best cost for first iter
    globalsolution=pop(indx,:); % Store best solution for first iter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start iterations
    iter = 0; %Initialization to the loop for the iterations
    while iter < rdpso.niter
       iter = iter + 1;
       alpha=alphamax-((alphamax-alphamin)/rdpso.niter)*iter;
       r1 = rand(rdpso.n,rdpso.ndim); % random numbers
       r2 = rand(rdpso.n,rdpso.ndim); % random numbers
       % The following to evaluate equation (15) in paper no. 1.
       phi= c1*r1./(c1*r1+c2*r2); 
       % The following to evaluate equation (17) in paper no. 1.
       Cmean=sum(localsolution,1); 
       Cmeant=(1/rdpso.n)*repmat(Cmean,rdpso.n,1);
       % The following to evaluate equation (14) in paper no. 1.
       localfocus=(phi.*localsolution)+(1-phi).*(repmat(globalsolution,rdpso.n,1));
    
       % Update the position and the velocity
       % The following to evaluate equation (19) in paper no. 1.
       vel = alpha*abs(Cmeant-pop).*normrnd(0,1,rdpso.n,rdpso.ndim)+1.5*(localfocus-pop);
       % The following to evaluate equation (13) in paper no. 1.
       pop = pop + vel;
       %pop(pop>100) = 100;
       %pop(pop<0) = 0;
       
       % Evaluate the total population at the current iteration.
       parfor i = 1:rdpso.n
        [cost(1,i)]=costfunction(pop(i,:),neuron,params,lambdatype,rdpso); % Results in number of pop cost values
         %a1 = cat(1,a1,a);
         %b1 = cat(1,b1,b);
       end
    
       % The following uses equations (10) and (11) in paper no. 1 to update the local or
       % personal best positions for the particles and the global best
       % position.
       bettercost = cost < localcost;
       localcost = localcost.*not(bettercost) + cost.*bettercost;
       localsolution(bettercost,:) = pop(bettercost,:);
       rdpsosol.costplot(jn,iter) = min(cost);
       % Update the global solution per each iteration
       [temp, t] = min(cost);
       if temp<globalcost
       globalsolution=pop(t,:); indx=t; globalcost=temp;
       end
    end %End of the iterations loop
    
 % Handeling the final results
 % The results will be available in the workspace in the matrix results,
 % where results=[global cost per run,x1,x2,...,xn, time per run]
    time=toc;
    if rdpso.nruns==1
        results (1,1)=globalcost;
        results (1,2:(rdpso.ndim+1))=globalsolution;
        results (1,(rdpso.ndim+2))=time;
    else
        results (jn,1)=globalcost;
        results (jn,2:(rdpso.ndim+1))=globalsolution;
        results (jn,(rdpso.ndim+2))=time;
    end
    %toc
    %jn
    rdpsosol.globalsolution = globalsolution;
    rdpsosol.globalcost = globalcost;
end % End of no of runs loop
end

function [z] = costfunction(electrodestim,neuron,params,lambdatype,rdpso) % Cost
% z = solution to be minimzied
% ec = all electrode variables

if lambdatype == 1
    ec = (electrodestim./100).*rdpso.currentub; % Electrical stimulation values
    eo = zeros(length(electrodestim),1);
elseif lambdatype == 2
    ec = zeros(length(electrodestim),1);
    eo = (electrodestim./100).*rdpso.optoub; % Electrode optical stimulation value
else
    ec = (electrodestim(1:100)./100).*rdpso.currentub;
    eo = (electrodestim(101:200)./100).*rdpso.optoub; % Electrode optical stimulation value
end

% Electrical / Optical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(params.numneurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(params.numneurons,1); % Initialize current vector
Il_Soma_Neurons = zeros(params.numneurons,1); % Initialize luminous intensity

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i);
    Il_Soma_Neurons = Il_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons.
Il_Neurons = Il_Soma_Neurons; % Summation of luminous intensity directly from stimulus.

% Calculate Lambda Hat
[lambdahat] = lamdacombinedfun(neuron,params,Ie_Neurons,Il_Neurons,lambdatype);

% Threshold Determination:
Neuron_Activated = lambdahat > neuron.lambdamod+1;

% Ratio Calculation

a = sum(Neuron_Activated(neuron.motion.number)); % Number active motion-tuned
b = sum(Neuron_Activated(neuron.nonmotion.number)); % Number of active non-motion

if a >= length(neuron.motion.number)*.25 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a; % Ratio of what we don't want to activate over what we want to activate
else
    z = inf;
end

end

function y = costfunction2(x)
y = sum(x.^2);
end