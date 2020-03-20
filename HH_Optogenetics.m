%% Section III
% Hodgkin Huxley model including current and optogenetics
clearvars; clc;
options = odeset('RelTol',1e-6); % Tolerances are used to limit the local discretization error, this is attempting to minimize rounding errors.
h = 0.01; % Step size

% HH Model Initial Conditions
t = [1:h:1000];
y = zeros(9,1); % Initial conditions calculated by running code with all 0 and letting it run

% HH Model 
%y = [-66.5;0.31762;0.052927;0.59466]; % Initial conditions calculated by running code with all 0 and letting it run
[time,y_sim] = ode45(@ydiff,[min(t) max(t)],y,options); % Calls ODE45 function with y, time, and options inputs
y_sol = y_sim; % Stores Y solution data seperately
figure; plot(time,y_sol(:,1),'linewidth',2); ylabel('Voltage mV'); xlabel('Time in ms'); title('HH Model, Iapp = 10mA, Pulsewidth = 3ms');
figure; plot(time,y_sol(:,2:4),'linewidth',2); legend('n','m','h'); ylabel('Activation Parameters'); xlabel('Time in ms'); title('HH Model, Iapp = 10mA, Pulsewidth = 3ms');
disp(['#steps taken by ode45: ' num2str(length(time))])

%% Functions

function dydt = ydiff(t,y)

dydt = zeros(9,1); % y is vector with values for : voltage, m,n,h
% Where m,n,h are voltage dependent gating variables that change the ion
% conductance values, gna,gk,gl, depending on how active that ion is at
% that voltage.
V = y(1); %  V is stored in the first column of y
n = y(2); %  n is stored in the second column of y
m = y(3); %  m is stored in the third column of y
h = y(4); %  h is stored in the fourth column of y
O1 = y(5);
O2 = y(6);
C1 = y(7);
C2 = y(8);
p = y(9);

% Constants
Ena = 50; % Ion Reversal potential, Voltage
Ek = -77; % Ion Reversal potential, Voltage
El = -54.4; % leak channel Reversal potential, Voltage
gna = 120; % Ion maximum conductance, mS/cm^3
gk = 36; % Ion maximum conductance, mS/cm^3
gl = 0.3; % leak channel maximum conductance, mS/cm^3
Cm = 1; % Capacitence across the membrane, uF/cm^2
Echr2 = 0; % reversal potential for ChR2
gchcr2 = 0.4; % max conductance (scaling)
G(V) = [10.6408-14.6408*exp(-V/42.7671)]/V ;% voltage-dependent rectification function
gamma = 0.1; % ratio of conductances of O2/O1
e12d = 0.011; 
e12 = e12d+c1*ln(1+I/c2); % rate constant for O1?O2 transition, I – irradiance
e21d = 0.008;
e21 =  e21d+c1 *ln(1+I/c2); % I – irradiance
Gd2 = 0.05; % rate constant for O2?C2 transition

if 5 <= t && t <= 8 % Current input, x, Iapp
    x = 10;
else
    x = 0;
end

% 6 constant equations - Pg. 23 & 24 textbook V = v
% These voltage dependent constants compose the values for the three gating variables
an = 0.01*(V+55)/(1-exp(-(V+55)/10)); % Alpha n
bn = 0.125*exp(-(V+65)/80); % Beta n
am = 0.1*(V+40)/(1-exp(-(V+40)/10));
bm = 4*exp(-(V+65)/18);
ah = 0.07*exp(-(V+65)/20);
bh = 1/(1+exp(-(V+35)/10));
Iion = gna*(m^3)*h*(V-Ena)+gk*(n^4)*(V-Ek)+gl*(V-El);

k1 = epsilon1*Fp;
k2 = epsilon2*Fp;
F = sigmaret*I*lambda/(wloss*hc); % photon flux: number of photons per molecule per second
s0(theta) = 0.5(1 + tanh(120(theta-0.1)));
G(V) = [10.6408-14.6408*exp(-V/42.7671)]/V;
Gd1 = 0.075+0.043 *tanh((V+20)/-20); % rate constant for O1?C1 transition
Gr = 4.34587 * 105 * exp(-0.0211539274*V); % rate constant for C2?C1 transition
Ichcr2 = gchcr2 * G(V) *(O1 + gamma*O2)*(V-Echr2);

% Diff eq:
dydt(1) = -Iion + Ichcr2 +x; % Change in voltage this step: Cm dV/dt = -(Iion + Ichr2) + Stimulus
dydt(2) = an*(1-n)-bn*n; % dn/dt, change in n gating variable this step
dydt(3) = am*(1-m)-bm*m; % dm/dt, change in m gating variable this step
dydt(4) = ah*(1-h)-bh*h; % dh/dt, change in h gating variable this step
dydt(5) = k1*C1-(Gd1+e12)*O1+e21*O2;
dydt(6) = k2*C2-(Gd2+e21)*O2+e12*O1;
dydt(7) = Gr + Gd1 * O1-k*C1;
dydt(8) = Gd2*O2-(k2+Gr)*C2;
dydt(9) = (S0 - p) / tauChR2;
end