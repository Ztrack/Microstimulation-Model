clearvars; clc;
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)

load('InitialConditionsFull.mat')
load('solrb1.mat');

%% Calculations

for i = 1:params.NumMotifs
    DistanceMotif(i) = mean(nonzeros(sqrt((motif.center.x(i)-motif.center.x).^2 + (motif.center.y(i) - motif.center.y).^2)));
end
mean(DistanceMotif) % Average distance between motifs

for i = 1:params.numneurons
    DistanceNeurons(i) = min(nonzeros(sqrt((neuron.x(i)-neuron.x).^2 + (neuron.y(i) - neuron.y).^2)));
end
mean(DistanceNeurons) % Average distance between neurons in motifs

NeuronStimDist = NaN(1,params.numneurons);
for i = 1:params.numneurons
    NeuronStimDist(i) = sqrt((neuron.x(i)-electrode.x(ElectrodeNo)).^2 + (neuron.y(i) - electrode.y(ElectrodeNo)).^2);
end

PadX= [480 1440 2400 3360 4320 480 1440 2400 3360 4320 480 1440 2400 3360 4320];
PadY= [800 800 800 800 800 2400 2400 2400 2400 2400 4000 4000 4000 4000 4000];
DistancePads = zeros(1,params.numpads);
for i = 1:params.numpads
    DistancePads(i,:) = sqrt((PadX(i)-PadX).^2 + (PadY(i) - PadY).^2); % Distance of each pad to every other pad from the center
end

for i = 1:params.numneurons
    DistanceNeurons(i) = min(nonzeros(sqrt((neuron.x(i)-neuron.x).^2 + (neuron.y(i) - neuron.y).^2)));
end

stepsol.current.all = zeros(numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
stepsol.current.excitatory = stepsol.current.all; stepsol.current.inhibitory = stepsol.current.all; stepsol.current.motion = stepsol.current.all; stepsol.current.nonmotion = stepsol.current.all;
for i = 1:length(solrb.e.I0)
    stepsol.current.all(:,i) = sum(solrb.e1                          <=solrb.e.I0(i),2)*100/params.numneurons;
    stepsol.current.excitatory(:,i) = sum(solrb.e1(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    stepsol.current.inhibitory(:,i) = sum(solrb.e1(:,neuron.inhibitory)     <=solrb.e.I0(i),2)*100/length(neuron.inhibitory);
    stepsol.current.motion(:,i) = sum(solrb.e1(:,neuron.motion.number)    <=solrb.e.I0(i),2)*100/length(neuron.motion.number);
    stepsol.current.nonmotion(:,i) = sum(solrb.e1(:,neuron.nonmotion.number) <=solrb.e.I0(i),2)*100/length(neuron.nonmotion.number);
    
    stepsol.opto.all(:,i) = sum(solrb.o1                          <=solrb.o.I0(i),2)*100/params.numneurons;
    stepsol.opto.excitatory(:,i) = sum(solrb.o1(:,neuron.excitatory)     <=solrb.o.I0(i),2)*100/length(neuron.excitatory);
    stepsol.opto.inhibitory(:,i) = sum(solrb.o1(:,neuron.inhibitory)     <=solrb.o.I0(i),2)*100/length(neuron.inhibitory);
    stepsol.opto.motion(:,i) = sum(solrb.o1(:,neuron.motion.number)    <=solrb.o.I0(i),2)*100/length(neuron.motion.number);
    stepsol.opto.nonmotion(:,i) = sum(solrb.o1(:,neuron.nonmotion.number) <=solrb.o.I0(i),2)*100/length(neuron.nonmotion.number);
end

% Options for plotting
options = struct;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = 1;
%% Current Solution Plots
options.x_axis = solrb.e.I0;

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.current.all,options); title('ICMS Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');
ylim([0 100]);
hold off
legend('off');

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 0 255]./255; % Blue : Excitatory
plot_areaerrorbar(stepsol.current.excitatory,options); title('ICMS Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');
hold on
options.color_area = [255 0 0]./255; % Red : Inhibitory
plot_areaerrorbar(stepsol.current.inhibitory,options);
title('ICMSExcitatory / Inhibitory Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');
%xline(find(mean(stepsol.current.excitatory) >= 50,1),'color','blue'); xline(find(mean(stepsol.current.inhibitory) >= 50,1),'color','red'); 
legend('Excitatory','Inhibitory','50% Excitatory','50% Inhibitory'); 
ylim([0 100]);
disp((find(mean(stepsol.current.excitatory) >= 50,1)-(find(mean(stepsol.current.inhibitory) >= 50,1)))*50);
hold off

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : Motion
plot_areaerrorbar(stepsol.current.motion,options);
hold on
options.color_area = [255 165 0]./255; % Orange : Non-Motion
plot_areaerrorbar(stepsol.current.nonmotion,options);
title('ICMS Motion Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons'); 
%xline(find(mean(stepsol.current.motion) >= 50,1),'color',[0 128 0]./255); xline(find(mean(stepsol.current.nonmotion) >= 50,1),'color',[255 165 0]./255);
legend('Motion','Non-Motion','50% Motion','50% Non-Motion');
ylim([0 100]);
%xlim([0 currentstop/h]);
disp((find(mean(stepsol.current.motion) >= 50,1)-(find(mean(stepsol.current.nonmotion) >= 50,1)))*50)
hold off

% Plot of % neurons activated for every pad

for i = 1:params.numpads
    pad.params.numneurons(i) = sum(neuron.pad == i); % Sums number of neurons per pad
end

stepsol.current.pads = zeros(params.numpads,numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
for j = 1:params.numpads
for i = 1:length(solrb.e.I0)
    stepsol.current.pads(j,:,i) = sum(solrb.e1(:,neuron.pad == j)                          <=solrb.e.I0(i),2)*100/pad.params.numneurons(j);
end
end

%Colormap for colors
map = [0 0 0; 239 122 37; 81 182 106; 50 127 183; 152 85 159; 225 152 153;]./255; % Colormap. First number is 0, then 1 to 15.
c1 = .85;
c2 = .7;
map = [map ; map(2,1:3)*c1 ; map(3,1:3)*c1 ; map(4,1:3)*c1 ; map(5,1:3)*c1 ; map(6,1:3)*c1];
map = [map ; 1 1 1 ; map(3,1:3)*c2 ; map(4,1:3)*c2 ; map(5,1:3)*c2 ; map(6,1:3)*c2];

% Cumulative pad summation
pad.away(1).num = 8; Cumulative_Activation(1,:,:) = stepsol.current.pads(pad.away(1).num,:,:); % 0 pads away
pad.away(2).num = [7 9 3 13]; Cumulative_Activation(2,:,:) = mean(stepsol.current.pads(pad.away(2).num,:,:)); % 1 pads away
pad.away(3).num = [2 4 12 14]; Cumulative_Activation(3,:,:) = mean(stepsol.current.pads(pad.away(3).num,:,:)); % 1.5 pads away
pad.away(4).num = [6 10]; Cumulative_Activation(4,:,:) = mean(stepsol.current.pads(pad.away(4).num,:,:)); % 2 pads away
pad.away(5).num = [1 5 11 15]; Cumulative_Activation(5,:,:) = mean(stepsol.current.pads(pad.away(5).num,:,:)); % 2.5 pads away

%Plotting figure for cumulative pad summation
options.handle     = figure; set(gcf,'Position',[000 000 800 700]); j = 1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('ICMS Pad Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons'); 
legend('Middle Pad','1 Pad Away','1.5 Pads Away','2 Pads Away','2.5 Pads Away');
ylim([0 100]); 
%xlim([0 currentstop/h]);
hold off

%% Optogenetics Solution Plots

options.x_axis = solrb.o.I0;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.opto.all,options); title('Optogenetics Activation Through Center Electrode'); xlabel('Irradiance mW/mm^2'); ylabel('Percentage Activated Neurons');
ylim([0 100]);
legend('off');
%xlim([0 currentstop/h]);
hold off


options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 0 255]./255; % Blue : Excitatory
plot_areaerrorbar(stepsol.opto.excitatory,options); title('Optogenetics Neuron Activation Through Center Electrode'); xlabel('Irradiance mW/mm^2'); ylabel('Percentage Activated Neurons');
hold on
options.color_area = [255 0 0]./255; % Red : Inhibitory
plot_areaerrorbar(stepsol.opto.inhibitory,options);
title('Optogenetics Excitatory / Inhibitory Activation Through Center Electrode'); xlabel('Irradiance mW/mm^2'); ylabel('Percentage Activated Neurons');
%xline(find(mean(stepsol.current.excitatory) >= 50,1),'color','blue'); xline(find(mean(stepsol.current.inhibitory) >= 50,1),'color','red'); 
legend('Excitatory','Inhibitory','50% Excitatory','50% Inhibitory'); 
ylim([0 100]);
disp((find(mean(stepsol.opto.excitatory) >= 50,1)-(find(mean(stepsol.opto.inhibitory) >= 50,1)))*50);
hold off

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : Motion
plot_areaerrorbar(stepsol.opto.motion,options);
hold on
options.color_area = [255 165 0]./255; % Orange : Non-Motion
plot_areaerrorbar(stepsol.opto.nonmotion,options);
title('Optogenetics Motion Activation Through Center Electrode'); xlabel('Irradiance mW/mm^2'); ylabel('Percentage Activated Neurons'); 
%xline(find(mean(stepsol.current.motion) >= 50,1),'color',[0 128 0]./255); xline(find(mean(stepsol.current.nonmotion) >= 50,1),'color',[255 165 0]./255);
legend('Motion','Non-Motion','50% Motion','50% Non-Motion');
ylim([0 100]);
%xlim([0 currentstop/h]);
disp((find(mean(stepsol.opto.motion) >= 50,1)-(find(mean(stepsol.opto.nonmotion) >= 50,1)))*50)
hold off

% Plot of % neurons activated for every pad

for i = 1:params.numpads
    pad.params.numneurons(i) = sum(neuron.pad == i); % Sums number of neurons per pad
end

stepsol.opto.pads = zeros(params.numpads,numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of opto
for j = 1:params.numpads
for i = 1:length(solrb.e.I0)
    stepsol.opto.pads(j,:,i) = sum(solrb.e1(:,neuron.pad == j)                          <=solrb.e.I0(i),2)*100/pad.params.numneurons(j);
end
end

%Colormap for colors
map = [0 0 0; 239 122 37; 81 182 106; 50 127 183; 152 85 159; 225 152 153;]./255; % Colormap. First number is 0, then 1 to 15.
c1 = .85;
c2 = .7;
map = [map ; map(2,1:3)*c1 ; map(3,1:3)*c1 ; map(4,1:3)*c1 ; map(5,1:3)*c1 ; map(6,1:3)*c1];
map = [map ; 1 1 1 ; map(3,1:3)*c2 ; map(4,1:3)*c2 ; map(5,1:3)*c2 ; map(6,1:3)*c2];

% Cumulative pad summation
pad.away(1).num = 8; Cumulative_Activation(1,:,:) = stepsol.opto.pads(pad.away(1).num,:,:); % 0 pads away
pad.away(2).num = [7 9 3 13]; Cumulative_Activation(2,:,:) = mean(stepsol.opto.pads(pad.away(2).num,:,:)); % 1 pads away
pad.away(3).num = [2 4 12 14]; Cumulative_Activation(3,:,:) = mean(stepsol.opto.pads(pad.away(3).num,:,:)); % 1.5 pads away
pad.away(4).num = [6 10]; Cumulative_Activation(4,:,:) = mean(stepsol.opto.pads(pad.away(4).num,:,:)); % 2 pads away
pad.away(5).num = [1 5 11 15]; Cumulative_Activation(5,:,:) = mean(stepsol.opto.pads(pad.away(5).num,:,:)); % 2.5 pads away

%Plotting figure for cumulative pad summation
options.handle     = figure; set(gcf,'Position',[000 000 800 700]); j = 1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('Optogenetics Pad Activation Through Center Electrode'); xlabel('Irradiance mW/mm^2'); ylabel('Percentage Activated Neurons'); 
legend('Middle Pad','1 Pad Away','1.5 Pads Away','2 Pads Away','2.5 Pads Away');
ylim([0 100]); 
%xlim([0 currentstop/h]);
hold off

%% Current Vs Optogenetics 1

% Plot Options
options.legend = 0;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation Through Center Electrode');

% Plot ICMS
options.color_area = [0 128 0]./255; % Green
options.x_axis = solrb.e.I0; % X-axis
plot_areaerrorbar(stepsol.current.all,options);
ylim([0 100]);

plot(0,'color',[110 38 158]./255); hold on;
legend('location','northeastoutside');
legend('ICMS','Optogenetics'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = options.color_area;
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); ylabel('Percentage Activated Neurons');
hold on;

% Plot Opto
options.color_area = [110 38 158]./255;
options.x_axis = solrb.o.I0;
plot_areaerrorbar(stepsol.opto.all,options);
ylim([0 100]);

% Axis Settings 2
ax2.XColor = options.color_area;
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 

%% ICMS and opto effects on excitatory vs. inhibitory cells

% Plot Options
options.legend = 0;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation Through Center Electrode');

% Plot ICMS
options.x_axis = solrb.e.I0; % X-axis
options.color_area = [0 0 255]./255 ; % Blue : Excitatory
plot_areaerrorbar(stepsol.current.excitatory,options);
hold on
options.color_area = [255 0 0]./255; % Red : Inhibitory
plot_areaerrorbar(stepsol.current.inhibitory,options);
hold on
ylim([0 100]);

plot(0,'color',[0 128 0]./255); hold on;
plot(0,'color',[255 165 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS Excitatory','ICMS Excitatory','Opto Excitatory','Opto Inhibitory'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = [0 128 0]./255;
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); ylabel('Percentage Activated Neurons');
hold on;

% Plot Opto
options.x_axis = solrb.o.I0;
options.color_area = [0 128 0]./255; % Excitatory
plot_areaerrorbar(stepsol.opto.excitatory,options);
hold on
options.color_area = [255 165 0]./255; % Inhibitory
plot_areaerrorbar(stepsol.opto.inhibitory,options); 
ylim([0 100]);

% Axis Settings
ax2.XColor = [110 38 158]./255;
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 

hold off

%% ICMS and opto effects on Motion vs. Nonmotion cells

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 51 0]./255; % Green : Motion
plot_areaerrorbar(stepsol.current.motion,options);
hold on
options.color_area = [179 119 0]./255; % Orange : Non-Motion
plot_areaerrorbar(stepsol.current.nonmotion,options);
hold on
options.color_area = [128 255 128]./255; % Green : Motion
plot_areaerrorbar(stepsol.opto.motion,options);
hold on
options.color_area = [255 217 179]./255; % Orange : Non-Motion
plot_areaerrorbar(stepsol.opto.nonmotion,options);
title('Motion Activation Through Center Electrode'); 
xlabel('Luminous Intensity mW/mm^2'); ylabel('Percentage Activated Neurons'); 
legend('MS Motion','MS Non-Motion','Opto Motion','Opto Non-Motion');
ylim([0 100]);
disp((find(mean(stepsol.current.motion) >= 50,1)-(find(mean(stepsol.current.nonmotion) >= 50,1)))*50)
hold off




%% Plots - Axon 'Towards' or 'Away' stimulation center
ElectrodeNo = 45;
neuron.toward = find(neuron.dirmult(:,ElectrodeNo) == 1); % Array of all 'toward' neurons
neuron.away = find(neuron.dirmult(:,ElectrodeNo) ~= 1);

stepsol.current.toward = zeros(numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
stepsol.current.away = stepsol.current.toward;
for i = 1:length(solrb.e.I0)
    stepsol.current.toward(:,i) = sum(solrb.e1(:,neuron.toward)     <=solrb.e.I0(i),2)*100/length(neuron.toward);
    stepsol.current.away(:,i) = sum(solrb.e1(:,neuron.away)     <=solrb.e.I0(i),2)*100/length(neuron.away);
end
options.x_axis = solrb.e.I0;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.current.toward,options); 
title('Neuron Activation Through Center Electrode'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');
ylim([0 100]);
hold on
options.color_area = [179 119 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.current.away,options);

legend('Toward','Away');

%% Plots Neuron RB Vs Distance

%Excitatory and Inhibitory RB vs distance
figure; set(gcf,'Position',[000 000 800 700]);
Excitatory_RB = Neuron_RB(1,neuron.excitatory);
Excitatory_Dist = NeuronStimDist(neuron.excitatory);
Inhibitory_RB = Neuron_RB(1,neuron.inhibitory);
Inhibitory_Dist = NeuronStimDist(neuron.inhibitory);
plot(Excitatory_Dist,Excitatory_RB,'.','color','blue');
hold on
plot(Inhibitory_Dist,Inhibitory_RB,'.','color','red'); 
hold on; 

x = 1:3500;
Y1 = 0.00644.*x.^2+6.771.*x--3580; % Model Fit Data
Y2 = 0.003956.*x.^2+0.5833.*x-68.35; % Model Fit Data
plot(Y1,'blue');
plot(Y2,'red');
%ylim([0 1.75e04]); xlim([0 3000]);
legend('Excitatory','Inhibitory'); 
title('Neuron RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Current (pA)');
hold off

% % Neuron RB Vs distance but with pad distinction
% figure; set(gcf,'Position',[000 000 800 700]); 
% for i = 1:params.numpads
%     a=[]; b=a; c=a;
%     a = find(neuron.pad == i);
%     b(1,:) = intersect(a,neuron.excitatory); b(2,:) = NeuronStimDist(b); b(3,:) = Neuron_RB(1,b(1,:)); b = sortrows(b',2);
%     c(1,:) = intersect(a,neuron.inhibitory); c(2,:) = NeuronStimDist(c); c(3,:) = Neuron_RB(1,c(1,:)); c = sortrows(c',2);
%     plot(b(:,2),b(:,3),'-','DisplayName',num2str(i));
%     hold on
%     plot(c(:,2),c(:,3),'--','DisplayName',num2str(i));
%     hold on
% end
% legend show; title('Neuron RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Current (pA)');
% ylim([0 1.75e04]); xlim([0 3000]);
% hold off; 

% Fraction of cells activated within and across pads that are activated from microstim in the mid pad (M3)
figure; set(gcf,'Position',[000 000 800 700]);
RB_Pad = zeros(params.numpads,1);
for i = 1:params.numpads
    RB_Pad(i) = nanmean(Neuron_RB(1,neuron.pad == i)); % Rheobase for each pad
    plot(DistancePads(8,i),RB_Pad(i),'.','DisplayName',num2str(i)); xlim([-50 3000]);
    hold on
end
legend show; title('Average Pad RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Current (pA)');
hold off; 



%% Bar plots for mean & PCT M3 Neurons

pct = [NaN,25,100];
for i = 1:5
cumneurons(i) = sum(pad.params.numneurons(pad.away(i).num));
end

for ii = 1:length(pct)
if ii == 1
    string = 'Current Stimulation = Average Chronaxi of M3 Neurons';
    CH_M3 = RB_Pad(8)*2;
else
    string = ['Current Stimulation = Chronaxi of ' num2str(pct(ii)) '% M3 Neurons'];
    CH_M3 = prctile(nanmean(Neuron_RB(:,neuron.pad == 8)),pct(ii))*2;
end

Neuron_Activated = (mean(Neuron_RB) < CH_M3); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',neuron.pad(Neuron_Activated)'];
b1 = a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13;
b2 = a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14;
b3 = a(:,2) == 6 | a(:,2) == 10;
b4 = a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15;
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);

% Plotting Number of neurons stimulated when x pads away
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away

y= [length(onepadsaway) length(onepointfivepadsaway) length(twopadsaway) length(twopointfivepadsaway)]; % Number Neurons activated
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated');
%y= [length(onepadsaway).*100/(cumneurons(2)) length(onepointfivepadsaway).*100/(cumneurons(3)) length(twopadsaway).*100/(cumneurons(4)*2) length(twopointfivepadsaway).*100/(cumneurons(5))]; % Number Neurons activated
%figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated');
y= [length(onepadsaway) length(onepointfivepadsaway) length(twopadsaway) length(twopointfivepadsaway)]/length(onepadsaway); % Number Neurons activated
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated');

% Plotting number of excitatory/inhibitory neurons stimulated when x
% pads away
y= [sum(neuron.type(onepadsaway) == 2) sum(neuron.type(onepadsaway) == 1);sum(neuron.type(onepointfivepadsaway) == 2) sum(neuron.type(onepointfivepadsaway) == 1);sum(neuron.type(twopadsaway) == 2) sum(neuron.type(twopadsaway) == 1); sum(neuron.type(twopointfivepadsaway) == 2) sum(neuron.type(twopointfivepadsaway) == 1)];
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Excitatory','Inhibitory');
y= [sum(neuron.type(onepadsaway) == 2) sum(neuron.type(onepadsaway) == 1);sum(neuron.type(onepointfivepadsaway) == 2) sum(neuron.type(onepointfivepadsaway) == 1);sum(neuron.type(twopadsaway) == 2) sum(neuron.type(twopadsaway) == 1); sum(neuron.type(twopointfivepadsaway) == 2) sum(neuron.type(twopointfivepadsaway) == 1)]/sum(neuron.type(onepadsaway) == 2);
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Excitatory','Inhibitory');

% Plotting number of motion / non motion tuned neurons stimulated when x
% pads away
Motion_Designation = ones(params.numneurons,1); Motion_Designation(neuron.motion.number) = 1; Motion_Designation(neuron.nonmotion.number) = 2;
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
%y= [sum(Motion_Designation(onepadsaway) == 1).*100/(cumneurons(2)*(.2)) sum(Motion_Designation(onepadsaway) == 2).*100/(cumneurons(2)*(.8));sum(Motion_Designation(onepointfivepadsaway) == 1).*100/(cumneurons(3)*(.2)) sum(Motion_Designation(onepointfivepadsaway) == 2).*100/(cumneurons(3)*(.8));sum(Motion_Designation(twopadsaway) == 1).*100/(cumneurons(4)*(.2)) sum(Motion_Designation(twopadsaway) == 2).*100/(cumneurons(4)*(.8));sum(Motion_Designation(twopointfivepadsaway) == 1).*100/(cumneurons(5)*(.2)) sum(Motion_Designation(twopointfivepadsaway) == 2).*100/(cumneurons(5)*(.8))];
%figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)]/sum(Motion_Designation(onepadsaway) == 1);
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');

end

%% Box plots for % Motion-Tuned

pct = [25];
for ii = 1:length(pct)

string = ['Current Stimulation = Chronaxi of ' num2str(pct(ii)) '% Motion Neurons'];
CH_M3 = prctile(nanmean(Neuron_RB(:,neuron.motion.number)),pct(ii))*2;

Neuron_Activated = (mean(Neuron_RB) < CH_M3); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',neuron.pad(Neuron_Activated)'];
b1 = a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13;
b2 = a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14;
b3 = a(:,2) == 6 | a(:,2) == 10;
b4 = a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15;
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);

% Plotting number of motion / non motion tuned neurons stimulated when x
% pads away
Motion_Designation = ones(params.numneurons,1); Motion_Designation(neuron.motion.number) = 1; Motion_Designation(neuron.nonmotion.number) = 2;
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
%y= [sum(Motion_Designation(onepadsaway) == 1).*100/(cumneurons(2)*(.2)) sum(Motion_Designation(onepadsaway) == 2).*100/(cumneurons(2)*(.8));sum(Motion_Designation(onepointfivepadsaway) == 1).*100/(cumneurons(3)*(.2)) sum(Motion_Designation(onepointfivepadsaway) == 2).*100/(cumneurons(3)*(.8));sum(Motion_Designation(twopadsaway) == 1).*100/(cumneurons(4)*(.2)) sum(Motion_Designation(twopadsaway) == 2).*100/(cumneurons(4)*(.8));sum(Motion_Designation(twopointfivepadsaway) == 1).*100/(cumneurons(5)*(.2)) sum(Motion_Designation(twopointfivepadsaway) == 2).*100/(cumneurons(5)*(.8))];
%figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)]/sum(Motion_Designation(onepadsaway) == 1);
figure; set(gcf,'Position',[000 000 800 700]); bar(x,y); title(string); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');

end

%% Distributions of synaptic weights

figure;
hist(neuron.connections(:,4)*100,50);
title(['Individual Synaptic Weight Contributions n= ' num2str(length(neuron.connections))]); xlabel('Synaptic Weight (% of FR)'); ylabel('Number of Hits'); 

synaptic_weight = zeros(1,length(neuron.type));
for i = 1:length(neuron.connections)
    synaptic_weight(neuron.connections(i,2)) = synaptic_weight(neuron.connections(i,2)) + neuron.connections(i,4)*100;

    if neuron.connections(i,3) == 1 % If bi-directional
        synaptic_weight(neuron.connections(i,1)) = synaptic_weight(neuron.connections(i,1)) + neuron.connections(i,4)*100;

    end
end
figure;
hist(synaptic_weight(1,neuron.excitatory),50);
title(['Summated Synaptic Weights Per Neuron n= ' num2str(length(neuron.excitatory))]); xlabel('Synaptic Weight (% of FR)'); ylabel('Number of Hits'); 

figure;
hist(synaptic_weight(1,neuron.motion.number),50);
title(['Summated Synaptic Weights Per Motion-Tuned Neuron n= ' num2str(length(neuron.motion.number))]); xlabel('Synaptic Weight (% of FR)'); ylabel('Number of Hits'); 

figure;
synaptic_weight_inhibitory = neuron.inhibitoryfactor*100;
hist(synaptic_weight_inhibitory,50);
title(['Synaptic Weight of Inhibitory Neurons n= ' num2str(length(neuron.inhibitory))]); xlabel('Synaptic Weight (% of FR)'); ylabel('Number of Hits'); 

%% Neuron RB Varying with Inhibitory effect multiplier
load('rbsola.mat');

for i = 1:params.NumMotifs
    a = find(neuron.motif == i);
    b = a(1);
    a = a(2:5);
    Motif_Excitatory_RB_A1(i) = mean(neuron.rb.a1(a));
    Motif_Excitatory_RB_A1_err(i) = var(neuron.rb.a1(a));
    Motif_Excitatory_RB_A2(i) = mean(neuron.rb.a2(a));
    Motif_Excitatory_RB_A2_err(i) = var(neuron.rb.a2(a));
    Motif_Excitatory_RB_A3(i) = mean(neuron.rb.a3(a));
    Motif_Excitatory_RB_A3_err(i) = var(neuron.rb.a3(a));
    Motif_Excitatory_RB_A4(i) = mean(neuron.rb.a4(a));
    Motif_Excitatory_RB_A4_err(i) = var(neuron.rb.a4(a));
    Motif_Excitatory_RB_A5(i) = mean(neuron.rb.a5(a));
    Motif_Excitatory_RB_A5_err(i) = var(neuron.rb.a5(a));
    Motif_Excitatory_RB_A6(i) = mean(neuron.rb.a6(a));
    Motif_Excitatory_RB_A6_err(i) = var(neuron.rb.a6(a));
    
    Motif_Inhibitory_RB(i) = mean(neuron.rb.a1(b));
    Motif_Inhibitory_RB_err(i) = var(neuron.rb.a1(b));
end

for i = 1:params.NumMotifs
    DistanceMotifStim(i) = mean(nonzeros(sqrt((motif.center.x(i)-electrode.y(ElectrodeNo)).^2 + (motif.center.y(i) - electrode.x(ElectrodeNo)).^2)));
end
% 
% figure; set(gcf,'Position',[000 000 800 700]); plot(DistanceMotifStim,Motif_Excitatory_RB_A1,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A2,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A3,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A4,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A5,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A6,'-');
% title('Motif Excitatory RB'); xlabel('Distance Motif To Stimulus'); ylabel('Applied Current AU'); legend('.01','.03','.1','0.125','0.15','.2');


% figure; set(gcf,'Position',[000 000 800 700]); plot(DistanceMotifStim,Motif_Excitatory_RB_A1-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A2-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A3-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A4-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A5-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A6-Motif_Inhibitory_RB,'-');
% title('Excitatory-Inhibitory RB Difference'); xlabel('Distance Motif To Stimulus'); ylabel('Applied Current AU'); legend('.01','.03','.1','0.125','0.15','.2');


% Options for plotting
options = struct;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

stepsol.current.all = zeros(numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
stepsol.current.excitatory = stepsol.current.all; stepsol.current.inhibitory = stepsol.current.all; stepsol.current.motion = stepsol.current.all; stepsol.current.nonmotion = stepsol.current.all; Neuron_Excited_Per_Step6 = stepsol.current.all;

for i = 1:length(solrb.e.I0)
    stepsol.current.all(:,i) = sum(neuron.rb.a1(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    stepsol.current.excitatory(:,i) = sum(neuron.rb.a2(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    stepsol.current.inhibitory(:,i) = sum(neuron.rb.a3(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    stepsol.current.motion(:,i) = sum(neuron.rb.a4(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    stepsol.current.nonmotion(:,i) = sum(neuron.rb.a5(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
    Neuron_Excited_Per_Step6(:,i) = sum(neuron.rb.a6(:,neuron.excitatory)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory);
end

options.color_area = [255 0 0]./255;
options.color_line = [255 0 0]./255;
plot_areaerrorbar(stepsol.current.all,options);
hold on
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(stepsol.current.excitatory,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(stepsol.current.inhibitory,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(stepsol.current.motion,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(stepsol.current.nonmotion,options);
hold on
options.color_area = [0 0 128]./255;
options.color_line = [0 0 128]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step6,options);

xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Neurons at Applied Currents'); xlabel('Current (pA)'); ylabel('Percentage Activated Neurons');
legend('delete','0.01','delete','0.03','delete','0.1','delete','0.125','delete','0.15','delete','0.2');


%% Neuron RB Varying with adding more Neurons
load('rbsolb.mat');

neuron.excitatory1 = find(repmat([0 1 1 1 1 1],1,params.NumMotifs) == 1);
neuron.excitatory2 = find(repmat([0 0 1 1 1 1 1],1,params.NumMotifs) == 1);
neuron.excitatory3 = find(repmat([0 0 0 1 1 1 1 1],1,params.NumMotifs) == 1);
neuron.excitatory4 = find(repmat([0 0 0 0 1 1 1 1 1],1,params.NumMotifs) == 1);
neuron.excitatory5 = find(repmat([0 0 0 0 0 1 1 1 1 1],1,params.NumMotifs) == 1);

% Options for plotting
options = struct;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

stepsol.current.all = zeros(numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
stepsol.current.excitatory = stepsol.current.all; stepsol.current.inhibitory = stepsol.current.all; stepsol.current.motion = stepsol.current.all; stepsol.current.nonmotion = stepsol.current.all; Neuron_Excited_Per_Step6 = stepsol.current.all; Neuron_Excited_Per_Step7 = stepsol.current.all; Neuron_Excited_Per_Step8 = stepsol.current.all;

for i = 1:length(solrb.e.I0)
    %stepsol.current.all(:,i) = sum(Neuron_RBB0                           <=solrb.e.I0(i),2)*100/length(Neuron_RBB0);
    stepsol.current.excitatory(:,i) = sum(neuron.rb.a1(:,neuron.excitatory1)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory1);
    stepsol.current.inhibitory(:,i) = sum(neuron.rb.a2(:,neuron.excitatory2)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory2);
    stepsol.current.motion(:,i) = sum(neuron.rb.a3(:,neuron.excitatory3)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory3);
    stepsol.current.nonmotion(:,i) = sum(neuron.rb.a4(:,neuron.excitatory4)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory4);
    Neuron_Excited_Per_Step6(:,i) = sum(neuron.rb.a5(:,neuron.excitatory5)     <=solrb.e.I0(i),2)*100/length(neuron.excitatory5);
end

options.color_area = [255 0 0]./255;
options.color_line = [255 0 0]./255;
%plot_areaerrorbar(stepsol.current.all,options);
hold on
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(stepsol.current.excitatory,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(stepsol.current.inhibitory,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(stepsol.current.motion,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(stepsol.current.nonmotion,options);
hold on
options.color_area = [0 0 128]./255;
options.color_line = [0 0 128]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step6,options);
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Excitatory Neurons at Applied Currents'); xlabel('Current (pA)'); ylabel('Percentage Activated Excitatory Neurons');
legend('delete','1 Inhibitory','delete','2 Inhibitory','delete','3 Inhibitory','delete','4 Inhibitory','delete','5 Inhibitory');

% Subtract 5-4-3-2-1 for 4 total plots
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(stepsol.current.excitatory-stepsol.current.inhibitory,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(stepsol.current.inhibitory-stepsol.current.motion,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(stepsol.current.motion-stepsol.current.nonmotion,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(stepsol.current.nonmotion-Neuron_Excited_Per_Step6,options);
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Excitatory Difference at Applied Currents'); xlabel('Current (pA)'); ylabel('Percentage Activated Excitatory Neurons');
legend('delete','1-2 Inhibitory','delete','2-3 Inhibitory','delete','3-4 Inhibitory','delete','4-5 Inhibitory','delete','5-6 Inhibitory');
ylim([0 100]);

%% Functions

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
        options.error      = 'c95';
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
    set(patch, 'HandleVisibility','off');
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    if options.legend == 1
    legend('location','northeastoutside');
    end

end