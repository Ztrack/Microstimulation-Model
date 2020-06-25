% Author: Haedong Kim
% Date: 2019-08-25
% It is an example of how to use HH_model function

%% run the model
% clear the environment
clear
close
clc

% set hyper-parameter I
I(1:100000) = 1000; 
%I(501:1000) = 0; 
step_size = 0.01;

% run the model
[V, g_Na, g_K] = HH_model(I, 0, step_size);
pks = findpeaks(V,'MinPeakProminence',50); % Finds the spikes
spikes = length(pks)
% plot the results
x_lb = sprintf('Time Unit %s sec', step_size);

% V voltage
figure(1)
plot(V)
legend('Voltage')
ylabel('mV')
xlabel(x_lb)
title('Action Potential over Time')
