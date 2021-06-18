%% Microstim Model Figure Plots

clearvars; clc;
load('InitialConditionsFull.mat');
load('solrb1.mat');
set(groot,'defaultLineLineWidth',4.0);
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesLineWidth',3);
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
options.error = 'std';
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
[~,motifsort] = sort(sqrt((params.sx/2-neuron.x).^2 + (params.sy/2-neuron.y).^2)); % Sort motifs by distance
for i = 1:params.NumMotifs
    if sum(intersect(neuron.motion.number,find(neuron.motif == motifsort(i)))) > 0 % If it has at least 1 motion-tuned neuron
        center_motif = motifsort(i); % Store as motif of interest
        break
    end
end

% Figure 1a - View of entire model
figure; set(gcf,'Position',[100 100 800 700]); map = [1 1 1]; set(gca,'FontWeight','bold');
imagesc(population.empty); colormap(map); hold on;
plot(x1,y1,'.','color','Black'); hold on;
plot(x2,y2,'.','color','red'); hold on; 
plot(x3,y3,'.','color','Blue'); hold on;
plot(x4,y4,'.','color','Green'); hold on;
axis square;
sizebox = 500; % sizebox from center to plot, in um to focus on
y = motif.center.x(center_motif);
x = motif.center.y(center_motif);
drawSquare(params,x,y,sizebox);
%title('Neural Population'); xlabel('X Position (mm)'); ylabel('Y Position (mm)'); 
legend('','Inhibitory Cell','Excitatory Cell Non-Motion Tuned','Excitatory Motion Tuned',''); legend('location','northwestoutside');
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);

% Figure 1b - Zoomed in view of model

figure; set(gcf,'Position',[100 100 800 700]); map = [1 1 1];
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
xlim([y-sizebox y+sizebox]);
ylim([x-sizebox x+sizebox]);
 
% Figure 1c - View of motif
sizebox = 100; % sizebox from center to plot, in um
drawSquare(params,x,y,sizebox);

figure; set(gcf,'Position',[100 100 800 700]); map = [1 1 1];
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
xlim([y-sizebox y+sizebox]);
ylim([x-sizebox x+sizebox]);

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
hleg = legend('ICMS','Optogenetics');
htitle = get(hleg,'Title'); set(htitle,'String','Stimulation')

% Axis Settings
axis square;
hold on
ax1 = gca; % current axes
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1_pos = ax1.Position; % position of first axes
xlabel('Current (mA)'); ylabel('% Activated Neurons');
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
subplot(2,1,1);
options.linespecification = '-';
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('ICMS Pad Activation'); xlabel('Current (mA)'); ylabel('% Activated Neurons'); 
legendx = string(round(unique(pad.distance(8,:))));
ylim([0 100]); 
%xlim([0 stopcurrent]);
hleg = legend('location','northeastoutside'); legend([append(legendx,'\mum')]); htitle = get(hleg,'Title'); set(htitle,'String','Pad Distance');
set(gca,'FontWeight','bold');
axis square;
%xlim([0 currentstop/h]);

% Compare pad distributions

% Figure 2B Part 2
subplot(2,1,2);
for i = 1:params.numpads
    pad.params.numneurons(i) = sum(neuron.pad == i); % Sums number of neurons per pad
end

stepsol.opto.pads = zeros(params.numpads,numrepeats,length(solrb.o.I0)); % Calculates the number of neurons excited at every step of opto
for j = 1:params.numpads
for i = 1:length(solrb.o.I0)
    stepsol.opto.pads(j,:,i) = sum(solrb.o1(:,neuron.pad == j)                          <=solrb.o.I0(i),2)*100/pad.params.numneurons(j);
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
%options.handle     = figure; set(gcf,'Position',[000 000 800 700]); 
options.x_axis = solrb.o.I0;
options.linespecification = '--';
j = 1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('Optogenetics Pad Activation'); xlabel('Luminous Intensity mW/mm^2'); ylabel('% Activated Neurons'); 

axis square;
hleg = legend('location','northeastoutside'); legend([append(legendx,'\mum')]); htitle = get(hleg,'Title'); set(htitle,'String','Pad Distance');
set(gca,'FontWeight','bold');
ylim([0 100]); 
xlim([0 stopopto]);
hold off

% Figure 2C - Remove Inhibitory neurons
% analysis: What happens if we remove inhibitory cells and what happens if we remove the inhibitory effect
% Plot of % neurons activated for every pad
options.x_axis = solrb.e.I0;
options.linespecification = '-';
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
title('Removed Inhibitory Neurons'); xlabel('Current (mA)'); ylabel('% Activated Neurons'); 
ylim([0 100]); 
%xlim([0 stopcurrent]);
hleg = legend('location','northeastoutside'); legend([append(legendx,'\mum')]); htitle = get(hleg,'Title'); set(htitle,'String','Pad Distance');
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
title('Removed Inhibition'); xlabel('Current (mA)'); ylabel('% Activated Neurons'); 
hleg = legend('location','northeastoutside'); legend([append(legendx,'\mum')]); htitle = get(hleg,'Title'); set(htitle,'String','Pad Distance');
ylim([0 100]); 
%xlim([0 stopcurrent]);
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
hleg = legend('location','northeastoutside'); legend('0-45\circ','','46-90\circ','','91-134\circ','','135-180\circ'); htitle = get(hleg,'Title'); set(htitle,'String','Angle to Electrode');
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
hleg = legend('location','northeastoutside'); legend('0-45\circ','','46-90\circ','','91-134\circ','','135-180\circ'); htitle = get(hleg,'Title'); set(htitle,'String','Angle to Electrode');
axis square;
set(gca,'FontWeight','bold');
Cumulative_Activation(1,:,:) = stepsol.opto.x1;
Cumulative_Activation(2,:,:) = stepsol.opto.x2;
Cumulative_Activation(3,:,:) = stepsol.opto.x3;
Cumulative_Activation(4,:,:) = stepsol.opto.x4;
[mean,slope] = comparedist(Cumulative_Activation)
%%
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

% Figure 4A1 - Plots Neuron RB Vs Distance
Neuron_RB = mean(solrb.e1);
figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,1);
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
Y1 = (-4.92).*X1 .^2 + (256.4).*X1 + -148.5; % Model Fit Data - excitatory
X2 = 0:max(Inhibitory_RB)/100:max(Inhibitory_RB);
Y2 = (-246.4).*X2 .^2 + (1892).*X2 + -382.5; % Model Fit Data - Inhibitory
plot(X2,Y2,'red'); hold on;
plot(X1,Y1,'blue'); hold on;
%xlim([0 stopcurrent]);
ylim([0 3500]);
title('ICMS Rheobase vs Distance'); ylabel('Distance to stimulus um'); xlabel('Current (mA)');
hleg = legend('location','northeastoutside'); legend('Inhibitory','Excitatory'); htitle = get(hleg,'Title'); set(htitle,'String','Neuron Type');
axis square;
set(gca,'FontWeight','bold');
hold off

% Figure 4A2 - Plots Neuron RB Vs Distance
Neuron_RB = mean(solrb.o1);
%figure; set(gcf,'Position',[000 000 800 700]);
subplot(2,1,2);
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
Y1 = (2.004e+07).*X1 .^2 +  (.882e+04).*X1 + -78.59; % Model Fit Data - excitatory
X2 = 0:max(Inhibitory_RB)/100:max(Inhibitory_RB);
%Y2 = (1.683e+06).*X2 .^ (1.187); % Model Fit Data - Inhibitory
Y2 = (-3.245e+07).*X2 .^2 + (9.946e+05).*X2 + -1191;
plot(X2,Y2,'red'); hold on;
plot(X1,Y1,'blue'); hold on;
%xlim([0 stopcurrent]);
ylim([0 3500]);
title('Optogenetic Rheobase vs Distance'); ylabel('Distance to stimulus um'); xlabel('Luminous Intensity mW/mm^2');
hleg = legend('location','northeastoutside'); legend('Inhibitory','Excitatory'); htitle = get(hleg,'Title'); set(htitle,'String','Neuron Type');
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
hleg = legend('location','northeastoutside'); legend('ICMS Excitatory','ICMS Inhibitory','Opto Excitatory','Opto Inhibitory');  htitle = get(hleg,'Title'); set(htitle,'String','Neuron Type');

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

% Axis Settings 2
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
hleg = legend('location','northeastoutside'); legend('ICMS','Optogenetics'); htitle = get(hleg,'Title'); set(htitle,'String','Stimulation');

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

hleg = legend('location','northeastoutside'); legend('ICMS Motion','ICMS Non-Motion','Opto Motion','Opto Non-Motion'); htitle = get(hleg,'Title'); set(htitle,'String','Stimulation/Type');

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

% Figure 5B- Double Y-plot of difference between motion/nonmotion
% Plot Options
options.legend = 0;
options.linespecification = '-';
options.handle     = figure; set(gcf,'Position',[000 000 800 700]);
title('ICMS Vs Optogenetics: Activation');

% Plot Difference in ICMS
options.x_axis = solrb.e.I0; % X-axis
options.color_area = [0 128 0]./255 ; % Green
plot_areaerrorbar(abs(stepsol.current.motion-stepsol.current.nonmotion),options);
hold on
ylim([0 100]);
xlim([0 stopcurrent]);

plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
plot([1 2],[0 0],'--','color',[0 128 0]./255); hold on;
hleg = legend('location','northeastoutside'); legend('ICMS','Optogenetics'); htitle = get(hleg,'Title'); set(htitle,'String','Stimulation');

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

%% fitdist
%data.e1.mean = fitdist(mean(solrb.e1)','Normal');
%data.o1.mean = fitdist(mean(solrb.o1)','Normal');

% Compare pad distributions
clear x1 x2 x3 mu slope
for i = 1:size(Cumulative_Activation,1) % Finger pad #
    x1 = squeeze(Cumulative_Activation(i,:,:))';
    for j = 1:size(x1,2)
        [mu(i,j),slope(i,j)] = cdffit(x1(:,j));
    end
end
for i = 1:size(Cumulative_Activation,1) % Finger pad # dist 1
    for j = 1:size(Cumulative_Activation,1)
        if i == j
            x2(i,j) = nan;
            x3(i,j) = nan;
        else
            x2(i,j) = sum(mu(i,:)<mu(j,:))/length(mu(i,:)); % How many times is dist 1 greater than dist 2?
            x3(i,j) = sum(slope(i,:)<slope(j,:))/length(slope(i,:));
        end
    end
    
end


%% Motion vs Nonmotion cdf compare
x1 = mean(stepsol.opto.motion)';
x2 = mean(stepsol.opto.nonmotion)';
x3 = mean(stepsol.current.motion)';
x4 = mean(stepsol.current.nonmotion)';
[mudiff,slopediff] = distdiff(x1,x2);

%%
figure;
x1 = mean(stepsol.current.motion)';
p1 = fitdist(x1,'Normal')
x=1:150; A=3;
y=A*exp(-(x-p1.mu).^2/(2*p1.sigma^2))+rand(size(x))*0.3;
plot(x,y); hold on;

p1cdf = normcdf(x1,p1.mu,p1.sigma);
f = polyfit(1:length(p1cdf),p1cdf,1);
sigma = f(1);
mu = normcdf(median(x1),p1.mu,p1.sigma);
y=A*exp(-(x-mu).^2/(2*sigma^2))+rand(size(x))*0.3;
plot(x,y);

%%
clear x1 x2 x3 mu slope
for i = 1:size(Cumulative_Activation,1)
    x1 = mean(squeeze(Cumulative_Activation(i,:,:)))';
    [mu(i),slope(i)] = cdffit(x1);
end
x1 = mean(squeeze(Cumulative_Activation(i,:,:)))';
x2 = mean(squeeze(Cumulative_Activation(2,:,:)))';
x3 = mean(squeeze(Cumulative_Activation(3,:,:)))';
x4 = mean(squeeze(Cumulative_Activation(4,:,:)))';

x1=x1(x1>20);
x2=x2(x2>20);
x3=x3(x3>20);
x4=x4(x4>20);

%%

x1 = mean(stepsol.current.inhibitory);
x2 = mean(stepsol.current.excitatory);
x3 = mean(stepsol.opto.inhibitory);
x4 = mean(stepsol.opto.excitatory);
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

function [mu,slope] = cdffit(x1)
p1 = fitdist(x1,'Normal');
%p1cdf = normcdf(x1,p1.mu,p1.sigma);
%f = polyfit(1:length(p1cdf),p1cdf,1);
%slope = f(1);
%mu = normcdf(median(x1),p1.mu,p1.sigma);
mu = p1.mu;
slope = p1.sigma;
end

function [mudiff,sigmadiff] = distdiff(x1,x2)
p1 = fitdist(x1,'Normal');
p1cdf = normcdf(x1,p1.mu,p1.sigma);
f = polyfit(1:length(p1cdf),p1cdf,1);
slope1 = f(1);
mu1 = normcdf(median(x1),p1.mu,p1.sigma);

p2 = fitdist(x1,'Normal');
p2cdf = normcdf(x2,p2.mu,p2.sigma);
f = polyfit(1:length(p2cdf),p2cdf,1);
slope2 = f(1);
mu2 = normcdf(median(x2),p2.mu,p2.sigma);

mudiff = mu2-mu1;
sigmadiff = slope2-slope1;
end

function [x2,x3] = comparedist(Cumulative_Activation)
% Compare axon distributions

for i = 1:size(Cumulative_Activation,1) % Finger pad #
    x1 = squeeze(Cumulative_Activation(i,:,:))';
    for j = 1:size(x1,2)
        [mu(i,j),slope(i,j)] = cdffit(x1(:,j));
    end
end
for i = 1:size(Cumulative_Activation,1) % Finger pad # dist 1
    for j = 1:size(Cumulative_Activation,1)
        if i == j
            x2(i,j) = nan;
            x3(i,j) = nan;
        else
            x2(i,j) = sum(mu(i,:)<mu(j,:))/length(mu(i,:)); % How many times is dist 1 greater than dist 2?
            x3(i,j) = sum(slope(i,:)<slope(j,:))/length(slope(i,:));
        end
    end
end
end
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

function imageReturn = drawSquare(params,x,y,size)
    imageReturn = image;
    intensity = 255;
    X1 = x - size;
    Y1 = y - size;
    X2 = x + size;
    Y2 = y + size;

    population.box.map = zeros(params.sx,params.sy);
    population.box.map(X1,Y1:Y2) = 1;
    population.box.map(X2,Y1:Y2) = 1;
    population.box.map(X1:X2,Y1) = 1;
    population.box.map(X1:X2,Y2) = 1;
    population.box.indices = find(population.box.map > 0);
    [yplot,xplot] = ind2sub([params.sx params.sy],population.box.indices); % Motion
    plot(xplot,yplot,'.','color','magenta'); hold on;
end

function F= fitmethis(data,varargin)
% FITMETHIS finds best-fitting distribution
%	F= fitmethis(X) fits all distributions available in MATLAB's function 
%	MLE to data in vector X and returns them ordered by either LogLikelihood 
%	or Akaike. Either 20 continuous or 5 discrete distributions
%	are used based on user input or the type of data supplied (see below)
%  A non-parametric (kernel) fitting is added, for data with "difficult" 
%	distributions. See comment below.
%
%	The function returns a structure array with fields: 
%	name: name of the distribution (see HELP MLE for a list)
%	par:  vector of parameter estimates (1, 2 or 3 values)
%	ci:   matrix of confidence limits, one column per parameter
%	LL:   Log-Likelihood of the data
%	aic:  Akaike Information Criterion
%
%	Optional arguments can be specified as name/value pairs. Argument 
%	names are case insensitive and partial matches are allowed.
%
%	Name			Value
%	'dtype'	character string that specifies if the data are continuous 
% 				('cont') or discrete ('disc'). If missing, the function 
% 				decides the data are discrete if all values of X are 
% 				natural numbers.
%
%	'ntrials' Specifies number of trials for the binomial distribution. 
% 				NTRIALS must be either a scalar or a vector of the same 
% 				size as X. If missing the binomial is not fitted.
%
%	'figure'	Either 'on' (default), or 'off'. If 'on' a plot of the 
% 				data and the best fitting distribution is produced (scaled 
%				to match the data). Requires aditional function 'plotfitdist'.
% 
%	'pdata'	Either 'on' (default) or 'off'. If 'on' a histogram of the
% 				data will be plotted.
%
%	'pdist'	Number of distributions to plot, starting with the best-fitting.
%				Default is 1 and maximum is 4. If a preferred distribution is chosen, 
%				only this will	be plotted. 
% 
%	'alpha'	A value between 0 and 1 specifying a confidence level
%				for CI of 100*(1-alpha) (default is 0.05).
%
%	'criterion'	Criterion used to order the fits. It can be: 'LL' for
%				Log-Likelihood (default), or 'AIC' for Akaike.
% 
%	'output'	If set to 'off' supresses output to the command window.
%				Default 'on'
%
%	'pref'	Name of the preferred distribution to plot
%
%	If X contains negative values, only the Normal distribution is fitted.
%	Also, if X contains values > 1 the Beta distribution is not fitted. If X
%	contains 0 some distributions are not fitted.
%
%	Requires Statistics Toolbox
% 
%	Example 1:
%	x= gamrnd(5,0.5,1000,1);
%	F= fitmethis(x);
%
%   Name      Par1       Par2        Par3         LogL
%  gamma   4.947e+00   5.034e-01               -1.461e+03
%    gev  -1.126e-02   8.965e-01   1.982e+00   -1.463e+03
%	  ...
% 
%	Example 2 (small bump to highlight kernel):
% 	x= gamrnd(5,0.5,1000,1);
% 	y= normrnd(5,0.3,80,1);
% 	F= fitmethis([x; y],'pdist',4)
warning off
% Defaults & Storage
dtype= {};
ntrials= [];
fig= 'on';
alpha= 0.05;
criterion= 'LL';
output= 'on';
F= struct('name',{},'par',[],'ci',[],'LL',[],'aic',[]);
prefdist= [];
plotdata= 'on';
plotdist= 1;
% Arguments
for j= 1:2:length(varargin)
	string= lower(varargin{j});
	switch string(1:min(4,length(string)))
		case 'dtyp'
			dtype= varargin{j+1};
		case 'ntri'
			ntrials= varargin{j+1};
		case 'figu'
			fig= varargin{j+1};
		case 'alph'
			alpha= varargin{j+1};
		case 'crit'
			criterion= varargin{j+1};
		case 'outp'
			output= varargin{j+1};
		case 'pref'
			prefdist= varargin{j+1};
		case 'pdat'
			plotdata= varargin{j+1};
		case 'pdis'
			plotdist= varargin{j+1};
		otherwise
			error('Unknown argument name');
	end
end
% Distributions
Cdist= {'normal'; 'exponential'; 'gamma'; 'logistic'; ...
		  'tlocationscale';...
		  'uniform'; 'ev'; 'rayleigh'; 'gev'; 'beta'; ...
		  'nakagami'; 'rician'; 'inversegaussian'; 'birnbaumsaunders'; ...
		  'gp'; 'loglogistic'; 'lognormal'; 'weibull'};
mustbepos= 11;
Ddist= {'binomial'; 'nbin'; 'unid';'geometric';'poisson'};
% Try determine data type: Discrete or Continuous (avoid 0)
if isempty(dtype)
	if isempty(find(1- (data+1)./(fix(data)+1), 1)) 
		dtype= 'disc';
	else
		dtype= 'cont';
	end
end
% Fit all... 
switch dtype(1:4)
	% Continuous
	case 'cont'
	for j= 1:numel(Cdist)
		
		% If negative values, only fit normal
		if min(data) < 0
			[phat,pci]= mle(data,'distribution','normal','alpha',alpha);
			F(j).name= Cdist{j};
			F(j).par= phat;
			F(j).ci=  pci;
			pdfv= pdf('normal',data,F(j).par(1),F(j).par(2));
			F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			F(j).aic= 2*2- 2*F(j).LL;
			break
		% Check: if values > 1 for Beta, do nothing
		elseif strcmp('beta',Cdist{j}) && max(data) > 1
			F(j).name= 'beta';
			F(j).LL= -Inf;
			F(j).aic= Inf;
		% Check: if values > 0 for some distr. (they are sorted), do nothing
		elseif  j >= mustbepos && min(data) == 0
			F(j).name= Cdist{j};
			F(j).LL= -Inf;
			F(j).aic= Inf;
		% Any other case do the fit ...
		else
			try
				[phat,pci]= mle(data,'distribution',Cdist{j},'alpha',alpha);
				F(j).name= Cdist{j};
				F(j).par= phat;
				F(j).ci=  pci;
				if numel(F(j).par) == 1
					pdfv= pdf(F(j).name,data,F(j).par(1));
				elseif numel(F(j).par) == 2
					pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2));
				else
					pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2),F(j).par(3));
				end
				F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
				F(j).aic= 2*numel(F(j).par)- 2*F(j).LL;
			catch
				F(j).name= Cdist{j};
				F(j).par= NaN;
				F(j).LL= -Inf;
				F(j).aic= Inf;
				fprintf('%s%s\n',Cdist{j},' distribution not applicable')
			end
		end
	end
	% Add non-parametric (kernel) estimation with ksdensity. The degrees 
	% of freedom (required for AIC) can be estimated as the trace of the 
	% smoothing matrix. This matrix is not returned by ksdensity, so I 
	% re-create it here assuming a normal distribution. This procedure 
	% returns pdf estimates identical to those returned by ksdensity. 
	% They are not required for anything, are calculated only as a test.
	% In general, kernel estimators will have higher LogLikelihoods but, 
	% due to the high number of parameters, they will have much higher AICs.
	% Use with caution.
	[Kpdf,mesh,bandw]= ksdensity(data,data);
	F(j+1).name= 'Kernel';
	F(j+1).par= bandw;
	F(j+1).LL=  sum(log(Kpdf));
	N= numel(data);
	w= zeros(N);
	S= zeros(N);
	for k= 1:N
		S(k,:)= normalf(data,data(k),bandw);
		w(k,:)= normalf(mesh,data(k),bandw);
	end
	Kpdf2= sum(w)/N; %#ok<*NASGU>
	dof= trace(S);
	F(j+1).aic= 2*dof- 2*F(j+1).LL;
	% Discrete
	case 'disc'
	for j= 1:numel(Ddist)
		% Binomial needs number of trials
		if strcmp('binomial',Ddist{j}) 
			F(j).name= 'binomial';
			if isempty(ntrials) || (numel(ntrials) > 1 && numel(data) ~= numel(ntrials))
				F(j).LL= -Inf;
				F(j).aic= Inf;
			else
				[phat,pci]= mle(data,'ntrials',ntrials,'distribution','binomial','alpha',alpha);
				F(j).par= phat;
				F(j).ci=  pci;
				pdfv= pdf('bino',data,ntrials,F(j).par(1));
				F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
				F(j).ntrials= ntrials;
			end
		else
			[phat,pci]= mle(data,'distribution',Ddist{j},'alpha',alpha);
			F(j).name= Ddist{j};
			F(j).par= phat;
			F(j).ci=  pci;
			if numel(F(j).par) == 1
				pdfv= pdf(F(j).name,data,F(j).par(1));
			elseif numel(F(j).par) == 2
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2));
			else
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2),F(j).par(3));
			end
			F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			F(j).aic= 2*numel(F(j).par)- 2*F(j).LL;
		end
	end
end
% Order by criterion
switch criterion
	case 'LL'
		index= sortrows([(1:size(F,2))',[F.LL]'],-2);
	case 'AIC'
		index= sortrows([(1:size(F,2))',[F.aic]'],2);
end
F= F(index(:,1));
% Nice screen output
if strcmp('on',output)
	fprintf('\n\t\t\t\tName\t\tPar1\t\tPar2\t\tPar3\t\tLogL\t\tAIC\n')
	for j= 1:size(F,2)
		switch numel(F(j).par)
			case 1
				fprintf('%20s \t%10.3e \t\t\t\t\t\t\t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
			case 2
				fprintf('%20s \t%10.3e \t%10.3e \t\t\t\t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
			case 3
				fprintf('%20s \t%10.3e \t%10.3e \t%10.3e \t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
		end
	end
end
% Choose the preferred distr. for plotting (if specified)
if ~isempty(prefdist)
	best= find(strcmp(cellstr(prefdist),cellstr({F.name})));
	if isempty(best)
		if strcmp('on',output) 
			fprintf('\n *** Warning: the preferred distr. is not in the list of fitted \n');
		end
	else
		plotfitdist(data,F(best),dtype,'pdat',plotdata);
		if strcmp('on',output)
			fprintf('\n *** Warning: Plotting preferred distr. [%s] NOT BEST FIT \n',F(best).name);
		end
	end
end
% Plot data histogram & best fits from 1 to plotdist
if strcmp('on',fig) && isempty(prefdist)
	plotfitdist(data,F,dtype,'pdist',plotdist,'pdat',plotdata);
end
% End of fitmethis
end

function y= normalf(x,mu,sigma)
% Normal function values at x with mean= mu and std= sigma
y= (1/(sigma*sqrt(2*pi)))* exp(-((x-mu).^2)./(2*sigma^2));
end

function plotfitdist(data,F,dtype,varargin)
%PLOTFITDIST plots data and fitted distribution
% PLOTFITDIST(DATA,F,DTYPE) plots one or more probability density 
% functions (specified in structure F), and the histogram 
% of values in vector DATA scaled by total area using 'trapz'. 
% 
% Structure array F must contain the following fields:
%  'name' (char. string) is the name of the distribution
%	(for examples of accepted names see 'HELP PDF').
%  'par' a vector of parameters of the distribution. The size of 'par'
%  (1 to 3) must match the requirements of function PDF for 
%  the specified distribution.
%  'F' may contain more than 1 distribution specification. Up to 4
%  distributions will be plotted, as specified in optional argument 'pdist'.
%  
%  DTYPE is a character string, either 'cont' for continuous 
%	distributions or 'disc' for discrete ones. Continuous 
%	distributions are plotted as a line, while discrete ones are 
%	plotted as an histogram superimposed on the data histogram.
% 
% If 'name' is 'binomial', structure F must have an additional 
%  field called 'ntrials', which must be a scalar specifying the number 
%  of trials. If 'ntrials' is a vector (i.e. different number of 
%  trials for each value in DATA) the plot cannot be made due to 
%  mismatch in the number of simulated data and number of trials.
% 
% Optional arguments:
% 'lwidth' Width of plotted line (default 1.5).
% 'legend' Wether to plot a legend. 'on' (default) or 'off'. The legend 
%          will be the names of each distribution plotted.
% 'pdata'  Wether to plot the data. 'on' (default), or 'off'.
% 'pdist'  Number of distributions to plot (1 to 4). Default 1. If the 
%          size F is smaller than pdist, it will be ignored.
% 
% The function is designed to work with 'fitmethis'. However, it can be
% used independently if the input arguments are provided as specified.
% Defaults
linewidth= 1.5;
legnd= 'on';
pdat= 'on';
pdist= 1;
linecol= [.6 .8 1; 1 .6 0; .6 .8 0; .8 .6 1];
% Arguments
for j= 1:2:length(varargin)
	string= lower(varargin{j});
	switch string(1:4)
		case 'lwid'
			linewidth= varargin{j+1};
		case 'lege'
			legnd= varargin{j+1};
		case 'pdat'
			pdat= varargin{j+1};
		case 'pdis'
			pdist= varargin{j+1};
		otherwise
			error(['Unknown argument name: ',varargin{j}]);
	end
end
% Histogram for Discrete/Continuous
if strcmp(dtype,'cont')
	x = min(data):range(data)/50:max(data);
	[bincount,binpos] = hist(data,min(50,numel(data)/5));
else
	x = unique(data);
	[bincount,binpos] = hist(data,x);
end
% Plot 1 or more distributions
figure(gcf); 
hold on
% Plot data first (just once)
bincount= bincount/trapz(binpos,bincount); % scaled frequencies
if strcmp(pdat,'on')
	datap= bar(binpos,bincount,'FaceColor',[.85 .85 .85],'EdgeColor',[1 1 1],'BarWidth',1);
end
set(gca,'Layer','top'); % Avoids lower edge of bars covering the axis. 
% Plot 'pdist' best distributions (up to 4 or size of F)
n2plot= min([pdist,4,size({F.name},2)]);
for j= 1:n2plot
	distname= F(j).name;
	par=      F(j).par;
	% Extract ntrials if it is there
	if isfield(F,'ntrials')
		ntrials= F(j).ntrials;
	end
	% Calculate predicted
	switch numel(par)
		case 1
			if strcmp('binomial',distname)
				y = pdf('bino',x,ntrials(1),par(1));
			elseif strcmp('Kernel',distname)
				y = ksdensity(data,x,'Bandwidth',par(1));
			else
				y = pdf(distname,x,par(1));
			end
		case 2
			y = pdf(distname,x,par(1),par(2));
		case 3
			y = pdf(distname,x,par(1),par(2),par(3));
	end
	% Plot distrib. pdf
	if strcmp('cont',dtype)
		model= plot(x,y,'r','LineWidth',linewidth,'Color',linecol(j,:));
	else
		model= bar(x,y,'FaceColor',linecol(j,:),'EdgeColor','none','BarWidth',0.5);
	end
end
xlabel('Data'); ylabel('PDF');
% Plot legend with the appropriate number of distr. names
if strcmp(legnd,'on')
	legend(['data' {F(1:n2plot).name}]); legend('boxoff'); 
end
end