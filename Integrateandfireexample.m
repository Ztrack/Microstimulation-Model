clearvars; clc;

%% Pyramidal Neuron Example
V_reset = -80e-3; % reset membrane potential (Overshoot)
V_e = -72.36e-3; % Resting Membrane Potential
V_th = -41.40e-3; % Threshold Membrane Potential
Rm = 160.47e6; % resistance across membrane
tau_m = 19.73e-3; % membrane time constant


dt = 0.0002; %Time step
T = 0:dt:1; % 1 second simulation

Vm(1) = V_reset;
Im = 2.165e-10;

for t=1:length(T)-1
    if Vm(t) > V_th
        Vm(t+1) = V_reset;
    else
        Vm(t+1) = Vm(t) + dt * ( -(Vm(t) - V_e) + Im * Rm) / tau_m;
    end
end

figure;
plot(T,Vm*1e3,'b-');
xlabel('Time(s)');
ylabel('Voltage (mV)');


% Creating F-I Curve
Im = 0:.0001e-9:1e-9;
S = zeros(1,length(Im));

for i = 1:length(Im)
    Vm(1) = V_reset;
    
    for t=1:length(T)-1
        if Vm(t) > V_th
            Vm(t+1) = V_reset;
            S(i) = S(i) + 1;
        else
            Vm(t+1) = Vm(t) + dt * ( -(Vm(t) - V_e) + Im(i) * Rm) / tau_m;
        end
    end
end

figure;
plot(Im,S);
title('Frequency-Input Curve');
xlabel('Voltage (V)');
ylabel('Firing Rate (Hz)');
%legend('Neocortex Pyramidal Layer 2/3','Neocortex Basket');
% 20Hz = 2.165e-10
%% Basket Neuron Example

V_reset = -80e-3; % reset membrane potential (Overshoot)
V_e = -67.52e-3; % Resting Membrane Potential
V_th = -40.25e-3; % Threshold Membrane Potential
Rm = 152.7e6; % Resistance across membrane
tau_m = 9.57e-3; % membrane time constant


dt = 0.0002;
T = 0:dt:1; % 1 second simulation

Vm(1) = V_reset;
Im = .2e-9;

for t=1:length(T)-1
    if Vm(t) > V_th
        Vm(t+1) = V_reset;
    else
        Vm(t+1) = Vm(t) + dt * ( -(Vm(t) - V_e) + Im * Rm) / tau_m;
    end
end

% Creating F-I Curve
Im = 0:.0001e-9:1e-9;
S = zeros(1,length(Im));

for i = 1:length(Im)
    Vm(1) = V_reset;
    
    for t=1:length(T)-1
        if Vm(t) > V_th
            Vm(t+1) = V_reset;
            S(i) = S(i) + 1;
        else
            Vm(t+1) = Vm(t) + dt * ( -(Vm(t) - V_e) + Im(i) * Rm) / tau_m;
        end
    end
end
% 40Hz = 0.1991e-9