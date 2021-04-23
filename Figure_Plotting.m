%% Microstim Model Figure Plots

clearvars; clc;
load('InitialConditionsFull.mat');
load('solrb1.mat');
set(groot,'defaultLineLineWidth',4.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)
%solrb.o1(2,:) = solrb.o1(1,:)+solrb.o1(1,:)*.02; solrb.e1(2,:) = solrb.e1(1,:)+solrb.e1(1,:)*.02;
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
options.linespecification = '-';
options.x_axis = solrb.e.I0;
stopcurrent = solrb.e.I0(find(stepsol.current.all(1,:) > 99,1));
startcurrent = solrb.e.I0(find(stepsol.current.all(1,:) > 5,1));
stopopto = solrb.o.I0(find(stepsol.opto.all(1,:) > 99,1));
startopto = solrb.o.I0(find(stepsol.opto.all(1,:) > 5,1));

%% Figure 1 - Model architecture

[y1,x1] = ind2sub([params.sx params.sy],population.axon.indices); % Axons
[y2,x2] = ind2sub([params.sx params.sy],population.inhibitory.indices); % Inhibitory
[y3,x3] = ind2sub([params.sx params.sy],population.excitatory.indices); % Non-Motion Excitatory
[y4,x4] = ind2sub([params.sx params.sy],population.motion.indices); % Motion
population.empty = zeros(params.sx,params.sy);

% Figure 1a - View of entire model
figure; set(gcf,'Position',[100 100 920 600]); map = [1 1 1]; set(gca,'FontWeight','bold');
imagesc(population.empty); colormap(map); hold on;
plot(x1,y1,'.','color','Black'); hold on;
plot(x2,y2,'.','color','red'); hold on; 
plot(x3,y3,'.','color','Blue'); hold on;
plot(x4,y4,'.','color','Green'); hold on;
axis square;
%title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); 
legend('','Inhibitory','Excitatory','Excitatory Motion Tuned'); legend('location','northwestoutside');
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);

% Figure 1b - Zoomed in view of model
center_motif = neuron.motif(find( sqrt((params.sx/2-neuron.x).^2 + (params.sy/2-neuron.y).^2) == min(sqrt((params.sx/2-neuron.x).^2 + (params.sy/2-neuron.y).^2))));
size = 500; % Size from center to plot, in um

figure; set(gcf,'Position',[100 100 675 600]); map = [1 1 1];
imagesc(population.empty); colormap(map); hold on;
plot(x1,y1,'.','color','Black'); hold on;
plot(x2,y2,'.','color','red'); hold on; 
plot(x3,y3,'.','color','Blue'); hold on;
plot(x4,y4,'.','color','Green'); hold on;
%title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); 
%legend('','Inhibitory','Excitatory','Motion'); legend('location','northeastoutside');
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
axis square;
xlim([motif.center.x(center_motif)-size motif.center.x(center_motif)+size]);
ylim([motif.center.y(center_motif)-size motif.center.y(center_motif)+size]);

% Figure 1c - View of motif
center_motif = neuron.motif(find( sqrt((params.sx/2-neuron.x).^2 + (params.sy/2-neuron.y).^2) == min(sqrt((params.sx/2-neuron.x).^2 + (params.sy/2-neuron.y).^2))));
size = 50; % Size from center to plot, in um

figure; set(gcf,'Position',[100 100 675 600]); map = [1 1 1];
imagesc(population.empty); colormap(map); hold on;
plot(x1,y1,'.','color','Black'); hold on;
plot(x2,y2,'.','color','red'); hold on; 
plot(x3,y3,'.','color','Blue'); hold on;
plot(x4,y4,'.','color','Green'); hold on;
%title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); 
%legend('','Inhibitory','Excitatory','Motion'); legend('location','northeastoutside');
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
axis square;
xlim([motif.center.x(center_motif)-size motif.center.x(center_motif)+size]);
ylim([motif.center.y(center_motif)-size motif.center.y(center_motif)+size]);

%% Figure 2 - Optogenetic & ICMS Activation

% Figure 2a - Current Vs Optogenetics 1
% Plot Options
options.legend = 0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation Through Center Electrode');

% Plot ICMS
options.color_area = [0 128 0]./255; % Green
options.x_axis = solrb.e.I0; % X-axis
plot_areaerrorbar(stepsol.current.all,options);
ylim([0 100]);
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS','Optogenetics'); 

% Axis Settings
axis square;
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('Percentage Activated Neurons');
set(gca,'FontWeight','bold');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); %ylabel('Percentage Activated Neurons');
hold on;

% Plot Opto
options.color_area = [0 128 0]./255; % Green
options.x_axis = solrb.o.I0;
options.linespecification = '--';
plot_areaerrorbar(stepsol.opto.all,options);
ylim([0 100]);
xlim([0 stopopto]);

% Axis Settings 2
axis square;
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
set(gca,'FontWeight','bold');

% Figure 2B - cumulative pad summation
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

% Plot 
options.x_axis = solrb.e.I0;
options.handle     = figure; set(gcf,'Position',[000 000 800 700]); j = 1;
options.linespecification = '-';
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('ICMS Pad Activation By Electrode Distance'); xlabel('Current (mA)'); ylabel('Percentage Activated Neurons'); 
legendx = string(round(unique(pad.distance(8,:))));
legend([append(legendx,'\mum')]);
ylim([0 100]); 
%xlim([0 stopcurrent]);
legend('location','northeastoutside');
set(gca,'FontWeight','bold');
axis square;
%xlim([0 currentstop/h]);
hold off

% Figure 2C - Remove Inhibitory neurons
% analysis: What happens if we remove inhibitory cells and what happens if we remove the inhibitory effect
% Plot of % neurons activated for every pad

for i = 1:params.numpads
    pad.params.numneurons(i) = sum(neuron.pad == i & neuron.type == 2); % Sums number of neurons per pad
end

stepsol.current.pads = zeros(params.numpads,numrepeats,length(solrb.e.I0)); % Calculates the number of neurons excited at every step of current
for j = 1:params.numpads
for i = 1:length(solrb.e.I0)
    stepsol.current.pads(j,:,i) = sum(solrb.e1(:,neuron.pad == j & neuron.type == 2)                          <=solrb.e.I0(i),2)*100/pad.params.numneurons(j);
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
options.handle     = figure; set(gcf,'Position',[000 000 800 700]); 
subplot(2,1,1);
j = 1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('Removed Inhibitory Neurons'); xlabel('Current (mA)'); ylabel('Percentage Activated Neurons'); 
ylim([0 100]); 
%xlim([0 stopcurrent]);
legend('location','northeastoutside'); 
legend([append(legendx,'\mum')]);
axis square;
set(gca,'FontWeight','bold');
%xlim([0 currentstop/h]);

% Figure 2D - Remove Inhibitory effect
% analysis: What happens if we remove inhibitory cells and what happens if we remove the inhibitory effect
% Plot of % neurons activated for every pad

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
%options.handle     = figure; set(gcf,'Position',[000 000 800 700]); 
subplot(2,1,2);
j = 1;
options.linespecification = '-';
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('Removed Inhibition'); xlabel('Current (mA)'); ylabel('Percentage Activated Neurons'); 
legend([append(legendx,'\mum')]);
ylim([0 100]); 
%xlim([0 stopcurrent]);
legend('location','northeastoutside');
%xlim([0 currentstop/h]);
axis square;
set(gca,'FontWeight','bold');

%% Figure 3 Calculations

%Find Axon Hillock - Electrode angle differences
neuron.axonhillockangletoelectrode = zeros(params.numelectrodes,params.numneurons);
for i = 1:params.numneurons
    x1 = neuron.x(i);
    y1 = neuron.y(i);
    
        for ii = 1:params.numelectrodes
            x2 = electrode.x(ii);
            y2 = electrode.y(ii);
            % tan(theta) = (y2-y1)/(x2-x1) % Finding theta, the angle from axon hillock to electrode
            theta1 = atan2d(y2-y1,x2-x1); % tan = Opposite/adjacent
            theta1 = mod(theta1,360); % Convert from [0 -180] to [0 360]
            theta2 = mod(theta1-Axon_Direction_Theta(i),360); % Angle Difference between electrode and axon hillock in degrees
            neuron.axonhillockangletoelectrode(i,ii) = min(360-theta2, theta2); % Normalize
        end
end

% Now we must bin electrodes into 45 degree bins (center on angle, so first chunk is -22.5 to 22.5deg)
chunksize = 45; % Size of bin in degrees
ElectrodeNo = 45;
neuron.anglebins = cell(180/chunksize,1);
for i = 1:length(neuron.excitatory)
    binnumber = ceil(neuron.axonhillockangletoelectrode(neuron.excitatory(i),ElectrodeNo)/45); % Round up to nearest integer
    neuron.anglebins{binnumber} = [neuron.anglebins{binnumber} neuron.excitatory(i)];
    
end

% Now we need to pair into groups of equal size(180/chunksize) dependeding
% on distance
neuron.distancetoelectrode = sqrt((electrode.x(ElectrodeNo)-neuron.x).^2+(electrode.y(ElectrodeNo)-neuron.y).^2);
maxdistance = 5; % Distance in um
neuronset = [];
x1 = []; x2 = [];
for i = 1:length(neuron.anglebins{1})
    x1(i,1) = neuron.anglebins{1}(i); % Neuron Numbers in matrix
    x2(i,1) = neuron.distancetoelectrode(x1(i,1)); % Neuron distance (in um) in matrix
    for ii = 2:length(neuron.anglebins)
        diff = min(abs(x2(i,1)-neuron.distancetoelectrode(neuron.anglebins{ii})));
        x1(i,ii) = find((neuron.distancetoelectrode==diff+x2(i,1)) | (neuron.distancetoelectrode==-diff+x2(i,1)),1); % Find neuron number for this close neuron
        x2(i,ii) = neuron.distancetoelectrode(x1(i,ii)); % Store it's distance to electrode
    end
end

for i = 1:length(x1)
    if any(abs(x2(i,1) - x2(i,:)) > maxdistance) % If any are true, does not meet condition. Remove from matrix
        x1(i,:) = nan; % Set flag to remove
        x2(i,:) = nan;
    end
end
x1(any(isnan(x1), 2), :) = [];
x2(any(isnan(x2), 2), :) = [];

%% Figure 3 - Effects of ICMS and optogenetic stimulation on axons

% Figure 3A - Toward-Away ICMS
stepsol.current.x1 = zeros(numrepeats,length(solrb.e.I0)); 
for i = 1:length(solrb.e.I0)
    stepsol.current.x1(:,i) = sum(solrb.e1(:,x1(:,1))     <=solrb.e.I0(i),2)*100/length(x1(:,1));
    stepsol.current.x2(:,i) = sum(solrb.e1(:,x1(:,2))     <=solrb.e.I0(i),2)*100/length(x1(:,2));
    stepsol.current.x3(:,i) = sum(solrb.e1(:,x1(:,3))     <=solrb.e.I0(i),2)*100/length(x1(:,3));
    stepsol.current.x4(:,i) = sum(solrb.e1(:,x1(:,4))     <=solrb.e.I0(i),2)*100/length(x1(:,4));
end

options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,1);
options.linespecification = '-';
options.color_area = [0 128 0]./255; % Green : All Neurons
options.x_axis = solrb.e.I0;
plot_areaerrorbar(stepsol.current.x1,options); 
xline(options.x_axis(find(stepsol.current.x1(1,:) > 50,1)),'color',options.color_area); hold on; options.color_area =[0 0 255]./255; % Blue
plot_areaerrorbar(stepsol.current.x2,options); hold on; 
xline(options.x_axis(find(stepsol.current.x2(1,:) > 50,1)),'color',options.color_area); options.color_area =[255 0 0]./255; % Red
plot_areaerrorbar(stepsol.current.x3,options); hold on; 
xline(options.x_axis(find(stepsol.current.x3(1,:) > 50,1)),'color',options.color_area); options.color_area =[255 165 0]./255; % Orange
plot_areaerrorbar(stepsol.current.x4,options); hold on;
xline(options.x_axis(find(stepsol.current.x4(1,:) > 50,1)),'color',options.color_area); 
ylim([0 100]);
xlim([startcurrent stopcurrent]);
title('Effects of ICMS Stimulation on Axons'); xlabel('Current (mA)'); ylabel('% Activated Neurons');
legend('0-45','','46-90','','91-134','','135-180');
legend('location','northeastoutside');
axis square;
set(gca,'FontWeight','bold');

% Figure 3B - Toward-Away  OPTO
stepsol.opto.x1 = zeros(numrepeats,length(solrb.o.I0)); 
for i = 1:length(solrb.o.I0)
    stepsol.opto.x1(:,i) = sum(solrb.o1(:,x1(:,1))     <=solrb.o.I0(i),2)*100/length(x1(:,1));
    stepsol.opto.x2(:,i) = sum(solrb.o1(:,x1(:,2))     <=solrb.o.I0(i),2)*100/length(x1(:,2));
    stepsol.opto.x3(:,i) = sum(solrb.o1(:,x1(:,3))     <=solrb.o.I0(i),2)*100/length(x1(:,3));
    stepsol.opto.x4(:,i) = sum(solrb.o1(:,x1(:,4))     <=solrb.o.I0(i),2)*100/length(x1(:,4));
end
%options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,2);
options.x_axis = solrb.o.I0;
options.linespecification = '--';
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.opto.x1,options); hold on; 
xline(options.x_axis(find(stepsol.opto.x1(1,:) > 50,1)),'color',options.color_area); options.color_area =[0 0 255]./255; % Blue
plot_areaerrorbar(stepsol.opto.x2,options); hold on; 
xline(options.x_axis(find(stepsol.opto.x2(1,:) > 50,1)),'color',options.color_area); options.color_area =[255 0 0]./255; % Red
plot_areaerrorbar(stepsol.opto.x3,options); hold on; 
xline(options.x_axis(find(stepsol.opto.x3(1,:) > 50,1)),'color',options.color_area); options.color_area =[255 165 0]./255; % Orange
plot_areaerrorbar(stepsol.opto.x4,options); hold on;
xline(options.x_axis(find(stepsol.opto.x4(1,:) > 50,1)),'color',options.color_area); 
ylim([0 100]);
xlim([startopto stopopto]);
title('Effects of Optogenetic Stimulation on Axons'); xlabel('Luminous Intensity mW/mm^2'); ylabel('% Activated Neurons');
legend('0-45','','46-90','','91-134','','135-180');
legend('location','northeastoutside');
axis square;
set(gca,'FontWeight','bold');

% Figure 3C - Difference plot
stepsol.current.xdiff = abs(stepsol.current.x1-stepsol.current.x4);
options.x_axis = solrb.e.I0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,1);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.current.xdiff,options); hold on; options.color_area =[0 0 255]./255; % Blue
ylim([0 100]);
xlim([0 stopcurrent]);
title('% Difference In Out of Phase ICMS'); xlabel('Current (mA)'); ylabel('% Activated Neurons');
legend('off');
axis square;
set(gca,'FontWeight','bold');

% Figure 3D - Opto Difference plot
stepsol.opto.xdiff = abs(stepsol.opto.x1-stepsol.opto.x4);
options.x_axis = solrb.o.I0;
options.linespecification = '--';
%options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,2);
options.color_area = [0 128 0]./255; % Green : All Neurons
plot_areaerrorbar(stepsol.opto.xdiff,options); hold on; options.color_area =[0 0 255]./255; % Blue
ylim([0 100]);
xlim([0 stopopto]);
title('% Difference In Out of Phase Optogenetics'); xlabel('Luminous Intensity mW/mm^2'); ylabel('% Activated Neurons');
legend('off');
axis square;
set(gca,'FontWeight','bold');

%% Figure 4 - Effects of ICMS and optogenetic stimulation on inhibitory vs excitatory cell types
% Figure 4A - Plots Neuron RB Vs Distance

Neuron_RB = mean(solrb.o1);
%Neuron_RB(:,neuron.away) = nan;
%Excitatory and Inhibitory RB vs distance
figure; set(gcf,'Position',[000 000 800 700]);
Excitatory_RB = Neuron_RB(1,neuron.excitatory);
Excitatory_Dist = NeuronStimDist(neuron.excitatory);
Inhibitory_RB = Neuron_RB(1,neuron.inhibitory);
Inhibitory_Dist = NeuronStimDist(neuron.inhibitory);
plot(Inhibitory_RB,Inhibitory_Dist,'.','color','red'); 
hold on; 
plot(Excitatory_RB,Excitatory_Dist,'.','color','blue');
hold on
set(gca,'FontWeight','bold');

X1 = 0:max(Excitatory_RB)/100:max(Excitatory_RB);
Y1 = (1.84e+07).*X1 .^2 +  (6.474e+04).*X1 + -190; % Model Fit Data - excitatory
X2 = 0:max(Inhibitory_RB)/100:max(Inhibitory_RB);
Y2 = (1.683e+06).*X2 .^ (1.187); % Model Fit Data - Inhibitory
plot(X2,Y2,'red'); hold on;
plot(X1,Y1,'blue'); hold on;
%xlim([0 stopcurrent]);
ylim([0 3500]);
legend('Inhibitory','Excitatory'); 
title('Neuron RB vs Distance to Stimulus'); ylabel('Distance to stimulus um'); xlabel('Current (mA)');
legend('location','northeastoutside');
axis square;
set(gca,'FontWeight','bold');
hold off

% Figure 4B ICMS and opto effects on excitatory vs. inhibitory cells

% Plot Options
options.legend = 0;
options.linespecification = '-';
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
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 0 255]./255); hold on;
plot([1 2],[0 0],'--','color',[255 0 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS Excitatory','ICMS Inhibitory','Opto Excitatory','Opto Inhibitory'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('% Activated Neurons');
axis square;
set(gca,'FontWeight','bold');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); %ylabel('% Activated Neurons');
hold on;

% Plot Opto
options.linespecification = '--';
options.x_axis = solrb.o.I0;
options.color_area = [0 0 255]./255 ; % Blue : Excitatory
plot_areaerrorbar(stepsol.opto.excitatory,options);
hold on
options.color_area = [255 0 0]./255; % Red : Inhibitory
plot_areaerrorbar(stepsol.opto.inhibitory,options); 
ylim([0 100]);
xlim([0 stopopto]);

% Axis Settings
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 
axis square;
set(gca,'FontWeight','bold');
hold off

% Figure 4C - Double Y-plot of difference inhib/excitatory

% Plot Options
options.legend = 0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation Through Center Electrode');

% Plot Difference in ICMS
options.x_axis = solrb.e.I0; % X-axis
options.color_area = [0 0 255]./255 ; % Blue : Excitatory
plot_areaerrorbar(abs(stepsol.current.inhibitory-stepsol.current.excitatory),options);
hold on
ylim([0 100]);
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 0 255]./255); hold on;
plot([1 2],[0 0],'--','color',[255 0 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS','Opto'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('% Activated Neurons');
axis square;
set(gca,'FontWeight','bold');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); %ylabel('% Activated Neurons');
hold on;

% Plot Opto
options.linespecification = '--';
options.x_axis = solrb.o.I0;
options.color_area = [0 0 255]./255 ; % Blue : Excitatory
plot_areaerrorbar(abs(stepsol.opto.inhibitory-stepsol.opto.excitatory),options);
ylim([0 100]);
xlim([0 stopopto]);

% Axis Settings
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 
axis square;
set(gca,'FontWeight','bold');
hold off

%% Figure 5 - ICMS stimulation on motion vs non motion tuned populations
% Figure 5A - Plots Neuron RB Vs Distance (motion vs nonmotion)

Neuron_RB = mean(solrb.o1);
%Neuron_RB(:,neuron.away) = nan;
%Excitatory and Inhibitory RB vs distance
figure; set(gcf,'Position',[000 000 800 700]);
Motion_RB = Neuron_RB(1,neuron.motion.number);
Excitatory_Dist = NeuronStimDist(neuron.motion.number);
Nonmotion_RB = Neuron_RB(1,intersect(neuron.nonmotion.number,neuron.excitatory));
Inhibitory_Dist = NeuronStimDist(intersect(neuron.nonmotion.number,neuron.excitatory));
plot(Nonmotion_RB,Inhibitory_Dist,'.','color','red'); 
hold on; 
plot(Motion_RB,Excitatory_Dist,'.','color','blue');
hold on
set(gca,'FontWeight','bold');

X1 = 0:max(Motion_RB)/100:max(Motion_RB);
Y1 = (1.84e+07).*X1 .^2 +  (6.474e+04).*X1 + -190; % Model Fit Data - excitatory
X2 = 0:max(Nonmotion_RB)/100:max(Nonmotion_RB);
Y2 = (1.683e+06).*X2 .^ (1.187); % Model Fit Data - Inhibitory
%plot(X2,Y2,'red'); hold on;
%plot(X1,Y1,'blue'); hold on;
%xlim([0 stopcurrent]);
ylim([0 3500]);
legend('Motion','Nonmotion'); 
title('Neuron RB vs Distance to Stimulus'); ylabel('Distance to stimulus um'); xlabel('Current (mA)');
legend('location','northeastoutside');
axis square;
set(gca,'FontWeight','bold');
hold off

% Figure 5B
% Plot Options
options.legend = 0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Motion Activation Through Center Electrode');

% Plot ICMS
options.x_axis = solrb.e.I0; % X-axis
options.color_area = [0 128 0]./255 ; % Green
plot_areaerrorbar(stepsol.current.motion,options);
hold on
options.color_area = [255 165 0]./255; % Orange
plot_areaerrorbar(stepsol.current.nonmotion,options);
hold on
ylim([0 100]);
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
plot([1 2],[0 0],'--','color',[255 165 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS Motion','ICMS Non-Motion','Opto Motion','Opto Non-Motion'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('% Activated Neurons');
axis square;
set(gca,'FontWeight','bold');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); %ylabel('% Activated Neurons');
hold on;

% Plot Opto
options.linespecification = '--';
options.x_axis = solrb.o.I0;
options.color_area = [0 128 0]./255 ; % Green
plot_areaerrorbar(stepsol.opto.motion,options);
hold on
options.color_area = [255 165 0]./255; % Orange
plot_areaerrorbar(stepsol.opto.nonmotion,options); 
ylim([0 100]);
xlim([0 stopopto]);

% Axis Settings
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 
axis square;
set(gca,'FontWeight','bold');

%% Figure 5B- Double Y-plot of difference between motion/nonmotion
% Plot Options
options.legend = 0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation Through Center Electrode');

% Plot Difference in ICMS
options.x_axis = solrb.e.I0; % X-axis
options.color_area = [0 128 0]./255 ; % Green
plot_areaerrorbar(abs(stepsol.current.motion-stepsol.current.nonmotion),options);
hold on
ylim([0 100]);
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
legend('location','northeastoutside');
legend('ICMS','Opto'); 

% Axis Settings
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('% Activated Neurons');
axis square;
set(gca,'FontWeight','bold');

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
xlabel('Luminous Intensity mW/mm^2'); %ylabel('% Activated Neurons');
hold on;

% Plot Opto
options.linespecification = '--';
options.x_axis = solrb.o.I0;
options.color_area = [0 128 0]./255 ; % Green
plot_areaerrorbar(abs(stepsol.opto.motion-stepsol.opto.nonmotion),options);
ylim([0 100]);
xlim([0 stopopto]);

% Axis Settings
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt); 
axis square;
set(gca,'FontWeight','bold');
hold off

%% Figure 6B - Optimization
% 
% 6. Cloud figures
% 	- Optimization solutions 
% 	- Optogenetics alone
% 	- MS alone
% 	- Combination w/ different colors
% 		- Gray for optical
% 		- Black for electrical
% 	- We want to show the nonlinearities. What would happen if we superimpose the independent vs when we use them at them at the same time.


%% Supplementary 1: SA1. SA,RA,Mixed neuron figure update
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
        options.linespecification = '-';
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
    plot(options.x_axis, data_mean, options.linespecification, ...
        'color', options.color_line, ...
        'LineWidth', options.line_width);
    if options.legend == 1
    legend('location','northeastoutside');
    end

end