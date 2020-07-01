clearvars; clc;

% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3191342/
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766297/

% Parameters
options = odeset('RelTol',1e-6); % Tolerances are used to limit the local discretization error, this is attempting to minimize rounding errors.
t = [0 1000];
%y = [17.1991 0.8462 0.95622];
y = [0 0 0 0]; % Initial conditions calculated by running code with all 0 and letting it run
type = 1;
x = 1;
y = [-64.0285 0.088919 0.032579 0.78245];

setGlobalx(x,type)
[time,y_sim] = ode45(@ydiff,[min(t) max(t)],y,options); % Calls ODE45 function with y, time, and options inputs
y_sol = y_sim; % Stores Y solution data seperately
pks = findpeaks(y_sol(:,1),'MinPeakHeight',-30,'MinPeakDistance',1); % Finds the spikes

figure(1); plot(time,y_sol(:,1),'linewidth',2); ylabel('Voltage mV'); xlabel('Time in ms'); title('HH Model');
figure(2); plot(time,y_sol(:,2:end),'linewidth',2); legend('n','m','h'); ylabel('Activation Parameters'); xlabel('Time in ms'); title('HH Model');


%%
x = linspace(0,30,10);

% Solutions
for i = 1:2
    
    type = i;
    if i == 1
        y = [-64.0285 0.088919 0.032579 0.78245];
    elseif i ==2
        y = [-65.0011 0.007849 0.016029 0.01459];
    end
    for ii = 1:length(x)
        setGlobalx(x(ii),type)
        [time,y_sim] = ode45(@ydiff,[min(t) max(t)],y,options); % Calls ODE45 function with y, time, and options inputs
        y_sol = y_sim; % Stores Y solution data seperately
        pks = findpeaks(y_sol(:,1),'MinPeakHeight',-30,'MinPeakDistance',1); % Finds the spikes
        spikes(i,ii) = length(pks);
    end
end
% Plotting
%figure(1); plot(time,y_sol(:,1),'linewidth',2); ylabel('Voltage mV'); xlabel('Time in ms'); title('HH Model');
%figure(2); plot(time,y_sol(:,2:end),'linewidth',2); legend('n','m','h'); ylabel('Activation Parameters'); xlabel('Time in ms'); title('HH Model');
%disp(['#steps taken by ode45: ' num2str(length(time))])
y1 = spikes(1,:);
y2 = spikes(2,:);
x(1) = nan; y1(1) = nan; y2(1) = nan;
figure; plot(x,y1); hold on; plot(x,y2); title('Frequency-Input Curve of Neurons'); xlabel('Applied Current (mA)'); ylabel('Frequency (Hz)'); legend('Inhibitory','Excitatory');

y3 = 31.71*sqrt(x--0.4279); y3(find(x<1.0911)) = 0;
y4 = 30.18.*sqrt(x-1.841); y4(find(x<2.2823)) = 0;
figure; plot(x,y1); hold on; plot(x,y2);
hold on;
plot(x,y3); hold on; plot(x,y4); 
title('Fit Frequency-Input Curve of Neurons'); xlabel('Applied Current (mA)'); ylabel('Frequency (Hz)'); legend('Inhibitory','Excitatory','Inhibitory Fit','Excitatory Fit');
ylim([0 max(y1)]);
xlim([0 max(x)]);
%% Functions

function dydt = ydiff(t,y)

dydt = zeros(4,1); % y is vector with values for : voltage, m,n,h
V = y(1);
n = y(2);
m = y(3);
h = y(4);
% Where m,n,h are voltage dependent gating variables that change the ion
% conductance values, gna,gk,gl, depending on how active that ion is at
% that voltage.

[x,type] = getGlobalx;

if t >= 0 % Current input, x, Iapp
    xin = x;
else
    xin = 0;
end

cm = 1; % Membrane capacitance

if type == 1 % Inhibitory Neuron
% Parameters
Ena = 55; % Ion Reversal potential, Voltage
Ek = -90; % Ion Reversal potential, Voltage
El = -65; % leak channel Reversal potential, Voltage
gna = 35; % Ion maximum conductance, mS/cm^2
gk = 9; % Ion maximum conductance, mS/cm^2
gl = 0.1; % leak channel maximum conductance, mS/cm^3

% These voltage dependent constants compose the values for the three gating variables
am = -0.1*(V+35)/(-1+exp(-(V+35)/10)); % Alpha n
bm = 4*exp(-(V+60)/18); % Beta n
ah = 0.07*exp(-(V+58)/20);
bh = 1/(1+exp(-(V+28)/10));
an = -0.01*(V+34)/(-1+exp(-(V+34)/10));
bn = 0.125*exp(-(V+44)/80);

elseif type == 2 % Excitatory Neuron
    
Ena = 60; % Ion Reversal potential, Voltage
Ek = -90; % Ion Reversal potential, Voltage
El = -65; % leak channel Reversal potential, Voltage
gna = 30; % Ion maximum conductance, mS/cm^3
gk = 100; % Ion maximum conductance, mS/cm^3
gl = 0.1; % leak channel maximum conductance, mS/cm^3

% These voltage dependent constants compose the values for the three gating variables
am = 0.182*(V+35)/(1-exp(-(V+35)/9)); % Alpha m
bm = -0.124*(V-35)/(1-exp((V-35)/9)); % Beta m
ah = 0.024*(V+50)/(1-exp(-(V+50)/5));
bh = -0.0091*(V-75)/(1-exp((V-75)/5));
an = 0.2*(V-20)/(1-exp(-(V-20)/9));
bn = -0.002*(V-20)/(1-exp((V-20)/9));
    
end

% Calculating Synaptic input
% where gsyn is the synaptic conductance and Esyn is the reversal
% potential of the synapse. For Esyn greater than the resting potential Vrest the synapse is depolarizing, i.e., excitatory, otherwise it is
% hyperpolarizing, i.e., inhibitory

%Isyn = Iampa+Inmda+Igaba;

% Differential Equations
dydt(1) = (-gna*(m^3)*h*(V-Ena)-gk*(n^4)*(V-Ek)-gl*(V-El)+xin)/cm; % Change in voltage this step
dydt(2) = (an*(1-n)-bn*n); % dn/dt, change in n gating variable this step
dydt(3) = (am*(1-m)-bm*m); % dm/dt, change in m gating variable this step
dydt(4) = (ah*(1-h)-bh*h); % dh/dt, change in h gating variable this step

end

function setGlobalx(val1,val2)
global x
x = val1;

global type
type = val2;
end

function [out1,out2] = getGlobalx
global x
out1 = x;

global type
out2 = type;
end
