clc; close all; clearvars;

%% Data Analyze
load('Pt4.mat');
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)

Neuron_Targeted = Neuron_Target(1,:); % Neuron #'s that we are targeting
Neuron_Non_Targeted = setdiff(1:NumNeurons,Neuron_Targeted); % Neuron #'s we are not targeting

x = [1 2];
y = [sum(Neuron_Activated(Neuron_Motion(1,:))) sum(Neuron_Activated(Neuron_NonMotion)); 0 0]; % Number Neurons activated
figure; bar(x,y); title(['Neurons Stimulated By Chronaxi From ' num2str(length(Stim_Electrode)) ' Electrodes']); ylabel('Neurons Activated'); legend('Motion-Tuned','Non-Motion'); xlim([0.5 1.5]);
y = [sum(Neuron_Activated(Neuron_Motion(1,:)))/sum(Neuron_Activated(Neuron_Motion(1,:))) sum(Neuron_Activated(Neuron_NonMotion))/sum(Neuron_Activated(Neuron_Motion(1,:))); 0 0]; % Number Neurons activated
figure; bar(x,y); title(['Neurons Stimulated By Chronaxi From ' num2str(length(Stim_Electrode)) ' Electrodes']); ylabel('Neurons Activated Normalized'); legend('Motion-Tuned','Non-Motion'); xlim([0.5 1.5]);


%% Inhibitory Excitation Optogenetics
% Activates Inhibitory Neurons with optogentics to inhibit excitatory
% neurons in the same motif
Inhibitory_Inactivation = zeros(NumNeurons,1);
for i = 1:NumMotifs
    if length(intersect(find(Neuron_Motif == i),Neuron_Motion(1,:))) == 0
        Inhibitory_Inactivation(find(Neuron_Motif == i)) = 1;
    end
end
Inhibitory_Inactivation = find(Inhibitory_Inactivation == 1); % Converts to Neuron #

%% Optogenetics optimization function

% For every electrode, increase optogenetics delivered until a motion
% neuron is inactivated (while loop) and apply optogenetic effect just
% before that point (if it is worthwhile). 

Opto_Thresh = 1;
I0 = 0:.1:405;
Opto_Stim = zeros(length(ElectrodeX),1);

for ii = 1:length(ElectrodeX) % For all electrodes, one electrode at a time
    b = [];
    for jj = 1:length(I0) % For all current steps, one at a time
        Opto_Soma_Neurons = I0_Soma_Neurons(:,ii).*I0(jj); % DIstance * current
        Neuron_Inactivated = Opto_Soma_Neurons >= Opto_Thresh; % All the neurons inactivated
        a = sum(Neuron_Inactivated(Neuron_Motion(1,:))); % All motion inactivated
        b(jj) = sum(Neuron_Inactivated(Neuron_NonMotion)); % all non-motion inactivated
        if a ~= 0 % if motion inactivated is not 0, break
            break
        end
    end
    if b(jj) > 0 % If we are inactivated non-motion, store this result
        Opto_Stim(ii) = I0(find(b == max(b),1));
    else
        Opto_Stim(ii) = 0; % otherwise, current for this electrode = 0
    end
end

%% Optogentics Application

Neuron_Activated_Post_Opto = Neuron_Activated; % Start with what was originally activated
Opto_Soma_Neurons = zeros(NumNeurons,1);

% Need array 600x1 of what is still active. If optop_thresh is reached, set
% Neuron_Activated_Post_Opto = 0
for i = 1:length(ElectrodeX) % Summation of stimulus for every neuron component by every electrode & its corresponding stimulus
    Opto_Soma_Neurons = Opto_Soma_Neurons + I0_Soma_Neurons(:,i).*Opto_Stim(i);
end
Neuron_Inactivated = Opto_Soma_Neurons > Opto_Thresh;
Neuron_Inactivated = find(Neuron_Inactivated == 1); % Converts to neuron#
Neuron_Activated_Post_Opto(Neuron_Inactivated) = 0;
Neuron_Activated_Post_Opto(Inhibitory_Inactivation) = 0;


%% Plotting
x = [1 2];
y = [sum(Neuron_Activated(Neuron_Motion(1,:))) sum(Neuron_Activated(Neuron_NonMotion)); sum(Neuron_Activated_Post_Opto(Neuron_Motion(1,:))) sum(Neuron_Activated_Post_Opto(Neuron_NonMotion));]; % Number Neurons activated
figure; bar(x,y); title(['Neurons Stimulated Post Optogenetics']); ylabel('Neurons Activated'); legend('Motion-Tuned','Non-Motion');
y = [sum(Neuron_Activated(Neuron_Motion(1,:)))/sum(Neuron_Activated(Neuron_Motion(1,:))) sum(Neuron_Activated(Neuron_NonMotion))/sum(Neuron_Activated(Neuron_Motion(1,:))); sum(Neuron_Activated_Post_Opto(Neuron_Motion(1,:)))/sum(Neuron_Activated_Post_Opto(Neuron_Motion(1,:))) sum(Neuron_Activated_Post_Opto(Neuron_NonMotion))/sum(Neuron_Activated_Post_Opto(Neuron_Motion(1,:)))]; % Number Neurons activated
figure; bar(x,y); title(['Neurons Stimulated Post Optogenetics']); ylabel('Neurons Activated Normalized'); legend('Motion-Tuned','Non-Motion');