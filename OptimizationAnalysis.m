clc; clear; close all;
load('InitialConditionsFull.mat');

n = 0.5; % Percent activation to count as 'activated'
x = 25; % The current/opto step must activate this many motion tuned neurons to store the step
load('SEoutputc.mat');
threshold.c = thresholdfun(output,units,n,x);
load('SEoutputo.mat');
threshold.o = thresholdfun(output,units,n,x);
%% Microstimulation + Optogenetics Optomization

neuron.inhibitoryfactor = 0.01;
neuron.threshold.rs = 29.0; % Value determined for 95%, used in testing
neuron.threshold.fs = 53.0; % % Value determined for 95%, used in testing
neuron.lambda = zeros(NumNeurons,1);
neuron.lambda(neuron.type == 1) = 40; % neuron.lambda for inhibitory Neurons
neuron.lambda(neuron.type == 2) = 20; % neuron.lambda for Excitatory Neurons

% lambdatype:
lambdatype = 1; % 3 = MS + Optogenetics (all excitatory - affects all cells indiscriminately)


%% PSO

if lambdatype == 1
    
    problem.nVar = 100;
    problem.VarMax = [threshold.c]';
    
elseif lambdatype == 2
    
    problem.nVar = 100;
    problem.VarMax = [threshold.o]';
    
else
    
    problem.nVar = 200;
    problem.VarMax = [threshold.c threshold.o]';
    
end

problem.CostFunction = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype); % Cost
problem.VarMin =  zeros(problem.nVar,1);

% Parameters of PSO
params.MaxIt = 50;        % Maximum Number of Iterations
params.AdaptiveItMax = 20; % Max number of adaptive iterations
params.AdaptiveItThreshold = .01; % if the absolute value of it - (it-1) is less than this, we say there is little enough change to end the search
params.nPop = 100;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = .99;        % Damping Ratio of Inertia Coefficient
params.c1 = 5;              % Personal Acceleration Coefficient
params.c2 = 5;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

out = PSOFunction2(problem, params);

%% Matlab PSO options

if lambdatype == 1
    nvars = 100;
    lb = zeros(1,nvars);
    ub = [threshold.c].*.5;
    
elseif lambdatype == 2
    nvars = 100;
    lb = zeros(1,nvars);
    ub = [threshold.o]*.5;
    
else
    
    nvars = 200;
    lb = zeros(1,nvars);
    ub = [threshold.c.*.5 threshold.o.*.5];
end

fun = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype);
options = optimoptions('particleswarm','SwarmSize',2000,'UseParallel',true,'display','iter');

%% Matlab PSO Once matlab
%ub = [];
[x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);


%% Matlab PSO Iterative

options = optimoptions('particleswarm','SwarmSize',1000,'UseParallel',true,'display','off');
numrepeats = 3;
CostsCurve = nan(5,numrepeats);

for i = 1:5
    
    lambdatype = i;
    fun = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype);
    if lambdatype == 1
        nvars = 100;
        lb = zeros(1,nvars);
        ub = [threshold.c]*.1;
        
    elseif lambdatype == 2
        nvars = 100;
        lb = zeros(1,nvars);
        ub = [threshold.o];
        
    else
        
        nvars = 200;
        lb = zeros(1,nvars);
        ub = [threshold.c*.1 threshold.o];
    end
    %ub = [];
    for ii = 1:numrepeats
        
        [x,fval,~,~] = particleswarm(fun,nvars,lb,ub,options);
        
        CostsCurve(i,ii) = fval;
        sol.type(i).repeat(ii).x = x;
        if nanmin(CostsCurve(i,:)) == fval % IF its the minimum for this set
            sol.type(i).bestx = x; % Save the solution
        end
        
    end
    
end

%% Matlab PSO Iterative Across #'s

for i = 1:5:250
    x = i;
    load('SEoutputc.mat');
    threshold.c = thresholdfun(output,units,n,x);
    load('SEoutputo.mat');
    threshold.o = thresholdfun(output,units,n,x);
    for ii = 1:4
        numrepeats = 1;
        for ii = 1:4
            
            lambdatype = ii;
            fun = @(electrodestim) MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype);
            for iii = 1:numrepeats
                [x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);
                CostsCurve(ii,iii,:) = fval;
            end
            
        end
    end
end
        
%% Plotting
CostsCurve(CostsCurve == inf) = NaN;
CostsCurve(CostsCurve == 0) = NaN;
figure;
for i = 1:size(CostsCurve,1)
    x = i;
    y = squeeze(nanmean(CostsCurve(i,:,:),2));
    yerr = squeeze(1.96.*nanstd(CostsCurve(i,:,:))/sqrt(numrepeats));
    errorbar(x,y,yerr);
    hold on;
end
xlim([0 8]);
ylim([0 6]);
xlabel('Lambda Calculation Type');
ylabel('Non-Motion / Motion Neuron Ratio');
title('Optimization Performance');
legend('MS','Opto','MS + Opto (+all)','MS + Opto (-all)','MS + Opto (+all) + Opto (-all)');

%% Plotting electrode stim maps

stimmap = zeros(4800,4800);
type = 1;
for i = 1:length(electrode.x)
    electrodemap = zeros(4800,4800);
    electrodemap(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1;
    electrodemap = bwdist(electrodemap) + electrodemap;
    electrodemap = sol.type(type).bestx(i)./electrodemap;
    stimmap = stimmap + electrodemap;
    
end
map = [1 1 1];
figure; imagesc(stimmap); set(gca,'YDir','normal');
hold on;
[y1,x1] = ind2sub([sx sy],population.inhibitory.indices); % Inhibitory
[y2,x2] = ind2sub([sx sy],population.excitatory.indices); % Non-Motion Excitatory
[y3,x3] = ind2sub([sx sy],population.motion.indices); % Motion
[y4,x4] = ind2sub([sx sy],population.axon.indices); % Axons
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;

%% Simulated Annealing
% https://www.mathworks.com/help/gads/simulannealbnd.html
options = optimoptions('simulannealbnd','display','iter');
x0 = ones(200,1);
[x,fval,exitflag,output] = simulannealbnd(fun,x0,lb,ub,options)

%% Direct Search
% https://www.mathworks.com/help/gads/patternsearch.html
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(100,1);
ub = [];
nonlcon = [];
x0 = rand(1,100).*threshold.o;
ub1 = [threshold.o];
%x0 = ub1.*rand(1,200);


options = optimoptions('patternsearch','display','iter');
[x,fval,exitflag,output] = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

%% Plot Iteration vs performance

% Options for plotting
options = struct;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green
plot_areaerrorbar(xsol1curve',options); 
title('Stimulation Performance');
xlabel('Iteration Number');
ylabel('Non-Motion / Motion-Tuned Neuron Ratio');
legend('Microstimulation');
ylim([0 7]);
hold off

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
plot_areaerrorbar(xsol1curve',options);
hold on;
options.color_area = [255 0 0]./255; % Red
plot_areaerrorbar(xsol4curve',options);
ylim([0 7]);
hold off
title('Stimulation Performance');
xlabel('Iteration Number');
ylabel('Non-Motion / Motion-Tuned Neuron Ratio');
legend('Microstimulation','Microstimulation + Optogenetics');

%% Plot Stimap
stimmap = zeros(4800,4800); 
sol = xsol1;
for i = 1:length(electrode.x)
    electrodemap = zeros(4800,4800);
    electrodemap(electrode.y(i)-electrode.radii:electrode.y(i)+electrode.radii,electrode.x(i)-electrode.radii:electrode.x(i)+electrode.radii) = 1;
    electrodemap = bwdist(electrodemap) + electrodemap;
    electrodemap = sol(i)./electrodemap;
    stimmap = stimmap + electrodemap;
    
end

%grad=colorGradient([1 1 1],[0 0 1],128); map = [0 0 0; grad];
grad=colorGradient([0 0 0],[1 1 1],2); map = [grad];
figure; imagesc(stimmap); colormap(map); colorbar;
set(gca,'YDir','normal'); set(gcf,'Position',[000 000 800 700]);
hold on;
[y1,x1] = ind2sub([sx sy],population.inhibitory.indices); % Inhibitory
[y2,x2] = ind2sub([sx sy],population.excitatory.indices); % Non-Motion Excitatory
[y3,x3] = ind2sub([sx sy],population.motion.indices); % Motion
[y4,x4] = ind2sub([sx sy],population.axon.indices); % Axons
plot(x4,y4,'.','color','Black')
hold on
plot(x1,y1,'.','color','red'); hold on; 
plot(x2,y2,'.','color','Blue'); hold on;
plot(x3,y3,'.','color','Green'); hold on;
title('Stimulation Intensity Map');



%%
function plot_areaerrorbar(data, options)
    options.color_line = options.color_area;
    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'std';
    end
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = mean(data,1);
    data_std  = std(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    end
    
    % Plotting the result
    figure(options.handle);
    x_vector = [options.x_axis', fliplr(options.x_axis')];
    patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    
    if options.legendswitch == 1
    legend(options.legend);
    end
    hold off;
    
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

%% Functions
function z = MotionRatioCombined(electrodestim,NumNeurons,neuron,lambdatype) % Cost
% z = solution to be minimzied
% ec = all electrode variables

if lambdatype == 1
    ec = electrodestim(1:100); % Electrical stimulation values
    eo = zeros(100,1);
elseif lambdatype == 2
    ec = zeros(100,1);
    eo = electrodestim(1:100); % Electrode optical stimulation value
else
    ec = electrodestim(1:100); % Electrical stimulation values
    eo = electrodestim(101:200); % Electrode optical stimulation values
end

% Electrical / Optical Field Summation

ElectrodeNo = 1:length(ec);
Ie_Axon_Neurons = zeros(NumNeurons,1); % Initialize current vector
Ie_Soma_Neurons = zeros(NumNeurons,1);
Il_Soma_Neurons = zeros(NumNeurons,1); % Initialize luminous intensity

for i = 1:length(ec) % Summation of current for every neuron component by every electrode & its corresponding current
    Ie_Axon_Neurons = Ie_Axon_Neurons + neuron.io.axon(:,ElectrodeNo(i)).*ec(i);
    Ie_Soma_Neurons = Ie_Soma_Neurons + neuron.io.soma(:,ElectrodeNo(i)).*ec(i).*neuron.dirmult(:,i);
    Il_Soma_Neurons = Il_Soma_Neurons + neuron.oo.soma(:,ElectrodeNo(i)).*eo(i);
end

Ie_Neurons = Ie_Soma_Neurons + Ie_Axon_Neurons; % Summation of current directly from stimulus + backpropogated up by axons.
Il_Neurons = Il_Soma_Neurons; % Summation of luminous intensity directly from stimulus.

% Calculate Lambda Hat
[lambdahat] = lamdacombinedfun(neuron,Ie_Neurons,Il_Neurons,neuron.inhibitoryfactor,lambdatype);

% neuron activation
NumTrials = 100;
simulation = 1;
dt = 1/1000;
bpct = .05;

% for i = 1:NumNeurons
%      output(1,i) = Simple_PoissonGen2(lambdahat(i), dt, NumTrials,simulation,bpct); % Calculate Lambda
% end
% Neuron_Activated = output > neuron.lambda+1;

% Threshold Method:
Neuron_Activated = zeros(NumNeurons,1);
Neuron_Activated(neuron.excitatory) = lambdahat(neuron.excitatory) > neuron.threshold.rs;
Neuron_Activated(neuron.inhibitory) = lambdahat(neuron.inhibitory) > neuron.threshold.fs;

% Ratio Calculation

a = sum(Neuron_Activated(neuron.motion.number)); % Number active motion-tuned = activated - inactivated
b = sum(Neuron_Activated(neuron.nonmotion.number)); % Number of active non-motion = activated - inactivated

if a >= length(neuron.motion.number)*.25 % If activated motion-tuned is less than threshold to create percept, ratio is no good
    z = b/a; % Ratio of what we don't want to activate over what we want to activate
else
    z = inf;
end

end