clearvars; clc;
load_Initial_Conditions = 1;
set(groot,'defaultLineLineWidth',2.0)
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesLineWidth',3)


%% Load Data

load('InitialConditionsFull.mat')
load('A1.mat');
Neuron_RB = Neuron_RBA1;
h = 50; % Step size
I0 = h:h:20000;
numrepeats = 100;
ElectrodeNo = 45;


h = 50; % Step size
I0 = h:h:25000;
numrepeats = 100;

%% Calculations

for i = 1:NumMotifs
    DistanceMotif(i) = mean(nonzeros(sqrt((MotifX(i)-MotifX).^2 + (MotifY(i) - MotifY).^2)));
end
mean(DistanceMotif) % Average distance between motifs

for i = 1:NumNeurons
    DistanceNeurons(i) = min(nonzeros(sqrt((NeuronX(i)-NeuronX).^2 + (NeuronY(i) - NeuronY).^2)));
end
mean(DistanceNeurons) % Average distance between neurons in motifs

PadX= [480 1440 2400 3360 4320 480 1440 2400 3360 4320 480 1440 2400 3360 4320];
PadY= [800 800 800 800 800 2400 2400 2400 2400 2400 4000 4000 4000 4000 4000];
DistancePads = zeros(1,NumPads);
for i = 1:NumPads
    DistancePads(i,:) = sqrt((PadX(i)-PadX).^2 + (PadY(i) - PadY).^2); % Distance of each pad to every other pad from the center
end

for i = 1:NumNeurons
    DistanceNeurons(i) = min(nonzeros(sqrt((NeuronX(i)-NeuronX).^2 + (NeuronY(i) - NeuronY).^2)));
end
%% Plots Neuron RB Vs Current

% Options for plotting
options = struct;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

options.handle     = figure;
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
Neuron_Excited_Per_Step1 = zeros(numrepeats,length(I0)); % Calculates the number of neurons excited at every step of current
Neuron_Excited_Per_Step2 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step3 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step4 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step5 = Neuron_Excited_Per_Step1;
for i = 1:length(I0)
    Neuron_Excited_Per_Step1(:,i) = sum(Neuron_RB                          <=I0(i),2)*100/NumNeurons;
    Neuron_Excited_Per_Step2(:,i) = sum(Neuron_RB(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step3(:,i) = sum(Neuron_RB(:,Neuron_Inhibitory)     <=I0(i),2)*100/length(Neuron_Inhibitory);
    Neuron_Excited_Per_Step4(:,i) = sum(Neuron_RB(:,Neuron_Motion(1,:))    <=I0(i),2)*100/length(Neuron_Motion);
    Neuron_Excited_Per_Step5(:,i) = sum(Neuron_RB(:,Neuron_NonMotion(:,1)) <=I0(i),2)*100/length(Neuron_NonMotion);
end
plot_areaerrorbar(Neuron_Excited_Per_Step1,options); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons');

options.handle     = figure;
options.color_area = [0 0 255]./255; % Blue : Excitatory
options.color_line = [0 0 255]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step2,options); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons');
hold on
options.color_area = [255 0 0]./255; % Red : Inhibitory
options.color_line = [255 0 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step3,options); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); 
title('Activated Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons');
xline(find(mean(Neuron_Excited_Per_Step2) >= 50,1),'color','blue'); xline(find(mean(Neuron_Excited_Per_Step3) >= 50,1),'color','red'); 
legend('','Excitatory','','Inhibitory','50% Inhibitory','50% Excitatory'); 
ylim([0 100]);
disp((find(mean(Neuron_Excited_Per_Step2) >= 50,1)-(find(mean(Neuron_Excited_Per_Step3) >= 50,1)))*50)
hold off

options.handle     = figure;
options.color_area = [0 128 0]./255; % Green : Motion
options.color_line = [0 128 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step4,options);
hold on
options.color_area = [255 165 0]./255; % Orange : Non-Motion
options.color_line = [255 165 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step5,options);
title('Activated Neurons at Applied Currents'); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons'); 
xline(find(mean(Neuron_Excited_Per_Step4) >= 50,1),'color',[0 128 0]./255); xline(find(mean(Neuron_Excited_Per_Step5) >= 50,1),'color',[255 165 0]./255);
legend('','Motion','','Non-Motion','50% Motion','50% Non-Motion');
ylim([0 100]);
disp((find(mean(Neuron_Excited_Per_Step4) >= 50,1)-(find(mean(Neuron_Excited_Per_Step5) >= 50,1)))*50)
hold off

% Plot of % neurons activated for every pad
Neuron_Pad_Des = zeros(NumPads,NumMotifsPerPad*NumNeuronsMotif);
for i = 1:NumPads
Neuron_Pad_Des(i,:) = find(Pad_Neuron == i);
end
Neuron_Excited_Per_Step_Per_Pad = zeros(NumPads,numrepeats,length(I0)); % Calculates the number of neurons excited at every step of current
for j = 1:NumPads
for i = 1:length(I0)
    Neuron_Excited_Per_Step_Per_Pad(j,:,i) = sum(Neuron_RB(:,Neuron_Pad_Des(j,:))                          <=I0(i),2)*100/40;
end
end

pad11 = squeeze(mean(Neuron_Excited_Per_Step_Per_Pad(11,:,:)));
map = [0 0 0; 239 122 37; 81 182 106; 50 127 183; 152 85 159; 225 152 153;]./255; % Colormap. First number is 0, then 1 to 15.
c1 = .85;
c2 = .7;
map = [map ; map(2,1:3)*c1 ; map(3,1:3)*c1 ; map(4,1:3)*c1 ; map(5,1:3)*c1 ; map(6,1:3)*c1];
map = [map ; 1 1 1 ; map(3,1:3)*c2 ; map(4,1:3)*c2 ; map(5,1:3)*c2 ; map(6,1:3)*c2];


n = 8; Cumulative_Activation(1,:,:) = Neuron_Excited_Per_Step_Per_Pad(n,:,:); % 0 pads away
n = [7 9 3 13]; Cumulative_Activation(2,:,:) = mean(Neuron_Excited_Per_Step_Per_Pad(n,:,:)); % 1 pads away
n = [2 4 12 14]; Cumulative_Activation(3,:,:) = mean(Neuron_Excited_Per_Step_Per_Pad(n,:,:)); % 1.5 pads away
n = [6 10]; Cumulative_Activation(4,:,:) = mean(Neuron_Excited_Per_Step_Per_Pad(n,:,:)); % 2 pads away
n = [1 5 11 15]; Cumulative_Activation(5,:,:) = mean(Neuron_Excited_Per_Step_Per_Pad(n,:,:)); % 2.5 pads away

options.handle     = figure;
j = 1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)*1.06),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
options.color_area = map(j+1,:); options.color_line = map(j+1,:); plot_areaerrorbar(squeeze(Cumulative_Activation(j,:,:)),options); hold on; j = j+1;
title('Activated Neurons at Different Pads'); xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons'); 
legend('','Middle Pad','','1 Pad Away','','1.5 Pads Away','','2 Pads Away','','2.5 Pads Away');
ylim([0 100]);
hold off

%% Plots Neuron RB Vs Distance

NeuronStimDist = NaN(1,NumNeurons);
for i = 1:NumNeurons
    NeuronStimDist(i) = sqrt((NeuronX(i)-ElectrodeX(ElectrodeNo)).^2 + (NeuronY(i) - ElectrodeY(ElectrodeNo)).^2);
end

% Excitatory and Inhibitory RB vs distance
% figure;
% Excitatory_RB = Neuron_RB(1,Neuron_Excitatory);
% Excitatory_Dist = NeuronStimDist(Neuron_Excitatory);
% Inhibitory_RB = Neuron_RB(1,Neuron_Inhibitory);
% Inhibitory_Dist = NeuronStimDist(Neuron_Inhibitory);
% plot(Excitatory_Dist,Excitatory_RB,'.','color','blue');
% hold on
% plot(Inhibitory_Dist,Inhibitory_RB,'.','color','red'); 
% hold on; load('fits.mat');
% plot(Excitatory_RB_Fit,'blue'); % Got from Cftools poly4
% plot(Inhibitory_RB_Fit,'red');
% ylim([0 1.75e04]); xlim([0 3000]);
% legend('Excitatory','Inhibitory','1398e-3*x^2 + 0.179*x + 15.74','299e-4*x^2 + 1.835*x + 587'); 
% title('Neuron RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Stimulus AU');
% hold off

% % Neuron RB Vs distance but with pad distinction
% figure; 
% for i = 1:NumPads
%     a=[]; b=a; c=a;
%     a = find(Pad_Neuron == i);
%     b(1,:) = intersect(a,Neuron_Excitatory); b(2,:) = NeuronStimDist(b); b(3,:) = Neuron_RB(1,b(1,:)); b = sortrows(b',2);
%     c(1,:) = intersect(a,Neuron_Inhibitory); c(2,:) = NeuronStimDist(c); c(3,:) = Neuron_RB(1,c(1,:)); c = sortrows(c',2);
%     plot(b(:,2),b(:,3),'-','DisplayName',num2str(i));
%     hold on
%     plot(c(:,2),c(:,3),'--','DisplayName',num2str(i));
%     hold on
% end
% legend show; title('Neuron RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Stimulus AU');
% ylim([0 1.75e04]); xlim([0 3000]);
% hold off; 

% Fraction of cells activated within and across pads that are activated from microstim in the mid pad (M3)
figure;
RB_Pad = zeros(NumPads,1);
for i = 1:NumPads
    RB_Pad(i) = nanmean(Neuron_RB(1,Pad_Neuron == i)); % Rheobase for each pad
    plot(DistancePads(8,i),RB_Pad(i),'.','DisplayName',num2str(i)); xlim([-50 3000]);
    hold on
end
legend show; title('Average Pad RB vs Distance to Stimulus'); xlabel('Distance to stimulus um'); ylabel('Stimulus AU');
hold off; 



%% Box plots 1
NumNeuronsPerPad = NumMotifsPerPad*NumNeuronsMotif;
CH_M3 = RB_Pad(8)*2; % Chronaxi for mid pad, M3
Neuron_Activated = (mean(Neuron_RB) < CH_M3); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',Pad_Neuron(Neuron_Activated)];
b1 = a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13;
b2 = a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14;
b3 = a(:,2) == 6 | a(:,2) == 10;
b4 = a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15;
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);


y= [length(onepadsaway) length(onepointfivepadsaway) length(twopadsaway) length(twopointfivepadsaway)]; % Number Neurons activated
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated');
y= [length(onepadsaway).*100/(NumNeuronsPerPad*4) length(onepointfivepadsaway).*100/(NumNeuronsPerPad*4) length(twopadsaway).*100/(NumNeuronsPerPad*2) length(twopointfivepadsaway).*100/(NumNeuronsPerPad*4)]; % Number Neurons activated
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away
figure; bar(x,y); title('Neurons Stimulated By Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated');

Neuron_Type2 = zeros(NumNeurons,1); Neuron_Type2(Neuron_Motion) = 1; Neuron_Type2(Neuron_NonMotion) = 2;
a = 100/(NumNeuronsPerPad*4*(4/5)); b = 100/(NumNeuronsPerPad*4*(1/5));
y= [sum(c(onepadsaway) == 1) sum(c(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 3);sum(c(onepointfivepadsaway) == 1) sum(c(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 3);sum(c(twopadsaway) == 1) sum(c(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 3);sum(c(twopointfivepadsaway) == 1) sum(c(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 3)];
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Excitatory Motion','Excitatory Non-Motion','Inhibitory','Inhibitory Oscillatory');
y= [sum(Neuron_Type(onepadsaway) == 2)*a sum(Neuron_Type2(onepadsaway) == 2)*b sum(Neuron_Type2(onepadsaway) == 3)*b;sum(Neuron_Type(onepointfivepadsaway) == 2)*a sum(Neuron_Type2(onepointfivepadsaway) == 2)*b sum(Neuron_Type2(onepointfivepadsaway) == 3)*b;sum(Neuron_Type(twopadsaway) == 2)*a sum(Neuron_Type2(twopadsaway) == 2)*b sum(Neuron_Type2(twopadsaway) == 3)*b;sum(Neuron_Type(twopointfivepadsaway) == 2)*a sum(Neuron_Type2(twopointfivepadsaway) == 2)*b sum(Neuron_Type2(twopointfivepadsaway) == 3)*b];
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Excitatory','Inhibitory','Inhibitory Oscillatory');

Motion_Designation = ones(NumNeurons,1); Motion_Designation(Neuron_Motion) = 1; Motion_Designation(Neuron_NonMotion) = 2;
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(onepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(onepointfivepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(onepointfivepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(twopadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(twopadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(twopointfivepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(twopointfivepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7))];
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)]/sum(Motion_Designation(onepadsaway) == 1);
figure; bar(x,y); title('Neurons Stimulated by Pad M3 Chronaxi'); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');

PadMatrix = zeros(sx,sy); px=sx/5; py=sy/3;
PadPercent = [y(4) y(2) y(1) y(2) y(4) y(3) y(1) 100 y(1) y(3) y(4) y(2) y(1) y(2) y(4)];
for i = 1:NumPads
    rangex = ((Padloc(i,3)-1)*px)+1:(Padloc(i,3)*px)-1;
    rangey = ((Padloc(i,2)-1)*py)+1:(Padloc(i,2)*py)-1;
    PadMatrix(rangey,rangex) = PadPercent(i);
end
figure; imagesc(PadMatrix); colorbar; title('Percent Neurons Activated By M3 Chronaxi'); xlabel('Position X'); ylabel('Position Y');

%% Box plots 2
% CH For 25% of M3 Neurons
CH_M3_25 = 325*2;
Neuron_Activated = (mean(Neuron_RB) < CH_M3_25); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',Pad_Neuron(Neuron_Activated)];
b1 = a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13;
b2 = a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14;
b3 = a(:,2) == 6 | a(:,2) == 10;
b4 = a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15;
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);

a = 100/(NumNeuronsPerPad*4); 
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away
y= [length(onepadsaway)*a length(onepointfivepadsaway)*a length(twopadsaway)*a length(twopointfivepadsaway)*a]; % Number Neurons activated
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated');
y= [length(onepadsaway) length(onepointfivepadsaway) length(twopadsaway) length(twopointfivepadsaway)]; % Number Neurons activated
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated');

a = 100/(NumNeuronsPerPad*4*(4/5)); b = 100/(NumNeuronsPerPad*4*(1/5));
y= [sum(Neuron_Type(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 3);sum(Neuron_Type(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 3);sum(Neuron_Type(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 3);sum(Neuron_Type(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 3)];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Excitatory','Inhibitory','Inhibitory Oscillatory');
y= [sum(Neuron_Type(onepadsaway) == 2)*a sum(Neuron_Type2(onepadsaway) == 2)*b sum(Neuron_Type2(onepadsaway) == 3)*b;sum(Neuron_Type(onepointfivepadsaway) == 2)*a sum(Neuron_Type2(onepointfivepadsaway) == 2)*b sum(Neuron_Type2(onepointfivepadsaway) == 3)*b;sum(Neuron_Type(twopadsaway) == 2)*a sum(Neuron_Type2(twopadsaway) == 2)*b sum(Neuron_Type2(twopadsaway) == 3)*b;sum(Neuron_Type(twopointfivepadsaway) == 2)*a sum(Neuron_Type2(twopointfivepadsaway) == 2)*b sum(Neuron_Type2(twopointfivepadsaway) == 3)*b];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Excitatory','Inhibitory','Inhibitory Oscillatory');

a = 100/(NumNeuronsPerPad*4*(.3)); b = 100/(NumNeuronsPerPad*4*(.7));
Motion_Designation = ones(NumNeurons,1); Motion_Designation(Neuron_Motion) = 1; Motion_Designation(Neuron_NonMotion) = 2;
y= [sum(Motion_Designation(onepadsaway) == 1)*a sum(Motion_Designation(onepadsaway) == 2)*b;sum(Motion_Designation(onepointfivepadsaway) == 1)*a sum(Motion_Designation(onepointfivepadsaway) == 2)*b;sum(Motion_Designation(twopadsaway) == 1)*a sum(Motion_Designation(twopadsaway) == 2)*b;sum(Motion_Designation(twopointfivepadsaway) == 1)*a sum(Motion_Designation(twopointfivepadsaway) == 2)*b];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)]/sum(Motion_Designation(onepadsaway) == 1);
figure; bar(x,y); title('Neurons Stimulated by CH of 25% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');

%% Box plots 3
% CH For 100% of M3 Neurons
CH_M3_100 = I0(25)*2;
Neuron_Activated = (mean(Neuron_RB) < CH_M3_100); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',Pad_Neuron(Neuron_Activated)];
b1 = a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13;
b2 = a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14;
b3 = a(:,2) == 6 | a(:,2) == 10;
b4 = a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15;
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);

y= [length(onepadsaway) length(onepointfivepadsaway) length(twopadsaway) length(twopointfivepadsaway)]; % Number Neurons activated
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated');
y= [length(onepadsaway).*100/(NumNeuronsPerPad*4) length(onepointfivepadsaway).*100/(NumNeuronsPerPad*4) length(twopadsaway).*100/(NumNeuronsPerPad*2) length(twopointfivepadsaway).*100/(NumNeuronsPerPad*4)]; % Number Neurons activated
x= [1 1.5 2 2.5]; % 1 - 2.5 pads away
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated');

c = zeros(NumNeurons,1); c(Neuron_Motion) = 1; c(Neuron_NonMotion) = 2;
a = 100/(NumNeuronsPerPad*4*(4/5)); b = 100/(NumNeuronsPerPad*4*(1/5));
y= [sum(c(onepadsaway) == 1) sum(c(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 2) sum(Neuron_Type2(onepadsaway) == 3);sum(c(onepointfivepadsaway) == 1) sum(c(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 2) sum(Neuron_Type2(onepointfivepadsaway) == 3);sum(c(twopadsaway) == 1) sum(c(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 2) sum(Neuron_Type2(twopadsaway) == 3);sum(c(twopointfivepadsaway) == 1) sum(c(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 2) sum(Neuron_Type2(twopointfivepadsaway) == 3)];
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Excitatory Motion','Excitatory Non-Motion','Inhibitory','Inhibitory Oscillatory');
y= [sum(Neuron_Type(onepadsaway) == 2)*a sum(Neuron_Type2(onepadsaway) == 2)*b sum(Neuron_Type2(onepadsaway) == 3)*b;sum(Neuron_Type(onepointfivepadsaway) == 2)*a sum(Neuron_Type2(onepointfivepadsaway) == 2)*b sum(Neuron_Type2(onepointfivepadsaway) == 3)*b;sum(Neuron_Type(twopadsaway) == 2)*a sum(Neuron_Type2(twopadsaway) == 2)*b sum(Neuron_Type2(twopadsaway) == 3)*b;sum(Neuron_Type(twopointfivepadsaway) == 2)*a sum(Neuron_Type2(twopointfivepadsaway) == 2)*b sum(Neuron_Type2(twopointfivepadsaway) == 3)*b];
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Excitatory','Inhibitory','Inhibitory Oscillatory');

Motion_Designation = ones(NumNeurons,1); Motion_Designation(Neuron_Motion) = 1; Motion_Designation(Neuron_NonMotion) = 2;
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(onepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(onepointfivepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(onepointfivepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(twopadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(twopadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7));sum(Motion_Designation(twopointfivepadsaway) == 1).*100/(NumNeuronsPerPad*4*(.3)) sum(Motion_Designation(twopointfivepadsaway) == 2).*100/(NumNeuronsPerPad*4*(.7))];
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)]/sum(Motion_Designation(onepadsaway) == 1);
figure; bar(x,y); title('Neurons Stimulated by CH of 100% M3 Neurons'); xlabel('Number of Pads Away'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');

%% Box plots 4
% CH For 25% of ALL Motion Neurons
CH_Motion_25 = prctile(nanmean(Neuron_RB(:,Neuron_Motion(1,:))),25)*2;
Neuron_Activated = (mean(Neuron_RB) < CH_Motion_25); % Stores if the neuron was activated for this applied stimulus
Neuron_Activated = find(Neuron_Activated == 1);
a = [Neuron_Activated',Pad_Neuron(Neuron_Activated)];
b1 = find(a(:,2) == 3 | a(:,2) == 7 | a(:,2) == 9 | a(:,2) == 13);
b2 = find(a(:,2) == 2 | a(:,2) == 4 | a(:,2) == 12 | a(:,2) == 14);
b3 = find(a(:,2) == 6 | a(:,2) == 10);
b4 = find(a(:,2) == 1 | a(:,2) == 5 | a(:,2) == 11 | a(:,2) == 15);
onepadsaway = a(b1);
onepointfivepadsaway = a(b2);
twopadsaway = a(b3);
twopointfivepadsaway = a(b4);


% a = 100/(NumNeuronsPerPad*4*(.3)); b = 100/(NumNeuronsPerPad*4*(.7));
% Motion_Designation = ones(NumNeurons,1); Motion_Designation(Neuron_Motion) = 1; Motion_Designation(Neuron_NonMotion) = 2;
% y= [sum(Motion_Designation(onepadsaway) == 1)*a sum(Motion_Designation(onepadsaway) == 2)*b;sum(Motion_Designation(onepointfivepadsaway) == 1)*a sum(Motion_Designation(onepointfivepadsaway) == 2)*b;sum(Motion_Designation(twopadsaway) == 1)*a sum(Motion_Designation(twopadsaway) == 2)*b;sum(Motion_Designation(twopointfivepadsaway) == 1)*a sum(Motion_Designation(twopointfivepadsaway) == 2)*b];
% figure; bar(x,y); title('Neurons Stimulated by CH of 25% ALL Motion Neurons'); xlabel('Number of Pads Away'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
% y= [sum(Motion_Designation(onepadsaway) == 1) sum(Motion_Designation(onepadsaway) == 2);sum(Motion_Designation(onepointfivepadsaway) == 1) sum(Motion_Designation(onepointfivepadsaway) == 2);sum(Motion_Designation(twopadsaway) == 1) sum(Motion_Designation(twopadsaway) == 2);sum(Motion_Designation(twopointfivepadsaway) == 1) sum(Motion_Designation(twopointfivepadsaway) == 2)];
% figure; bar(x,y); title('Neurons Stimulated by CH of 25% ALL Motion Neurons'); xlabel('Number of Pads Away'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');

a = 100/(length(Neuron_Motion)); b = 100/(NumNeurons-length(Neuron_Motion));
x = [1 2];
y= [length(intersect(Neuron_Activated,Neuron_Motion(1,:)))*a length(intersect(Neuron_Activated,Neuron_NonMotion(:,1)))*b; 0 0];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% ALL Motion Neurons'); xlabel('All Pads'); ylabel('Percent Neurons Stimulated'); legend('Motion','Non-Motion');
y= [length(intersect(Neuron_Activated,Neuron_Motion(1,:))) length(intersect(Neuron_Activated,Neuron_NonMotion(:,1))); 0 0];
figure; bar(x,y); title('Neurons Stimulated by CH of 25% ALL Motion Neurons'); xlabel('All Pads'); ylabel('Number Neurons Stimulated'); legend('Motion','Non-Motion');
y= [length(intersect(Neuron_Activated,Neuron_Motion(1,:))) length(intersect(Neuron_Activated,Neuron_NonMotion(:,1))); 0 0]/length(intersect(Neuron_Activated,Neuron_Motion(1,:)));
figure; bar(x,y); title('Neurons Stimulated by CH of 25% ALL Motion Neurons'); xlabel('All Pads'); ylabel('Neurons Stimulated Normalized'); legend('Motion','Non-Motion');


%% Temp 
% Plotting Ratio
% Equation = (# of non-motion tuned) / (# of motion tuned)

y= [90 279 ; 65 150];
figure; bar(x,y); title('Non-Optimized Vs Optimized Microstimulation'); xlabel('Non-Optimized Optimized'); ylabel('Number of Activated Neurons'); legend('Motion','Non-Motion');
xlim([0.5 2.5]);

y= [3.1 150/65 ; 0 0];
figure; bar(x,y); title('Non-Optimized Vs Optimized Microstimulation'); xlabel('Stimulation Strategy'); ylabel('Ratio Non-Motion To Motion'); legend('Optimized','Non-Optimized');
xlim([0.5 1.5]);

y= [150/65 65/45 ; 0 0];
figure; bar(x,y); title('Microstimulation Vs Microstimulation + Optogenetics'); xlabel('Stimulation Strategy'); ylabel('Ratio Non-Motion To Motion'); legend('MS Only','MS + Opto');
xlim([0.5 1.5]);




%% Neuron RB Varying with Inhibitory effect multiplier
load('A1.mat'); load('A2.mat'); load('A3.mat'); load('A4.mat'); load('A5.mat'); load('A6.mat');

for i = 1:NumMotifs
    a = find(Neuron_Motif == i);
    b = a(1);
    a = a(2:5);
    Motif_Excitatory_RB_A1(i) = mean(Neuron_RBA1(a));
    Motif_Excitatory_RB_A1_err(i) = var(Neuron_RBA1(a));
    Motif_Excitatory_RB_A2(i) = mean(Neuron_RBA2(a));
    Motif_Excitatory_RB_A2_err(i) = var(Neuron_RBA2(a));
    Motif_Excitatory_RB_A3(i) = mean(Neuron_RBA3(a));
    Motif_Excitatory_RB_A3_err(i) = var(Neuron_RBA3(a));
    Motif_Excitatory_RB_A4(i) = mean(Neuron_RBA4(a));
    Motif_Excitatory_RB_A4_err(i) = var(Neuron_RBA4(a));
    Motif_Excitatory_RB_A5(i) = mean(Neuron_RBA5(a));
    Motif_Excitatory_RB_A5_err(i) = var(Neuron_RBA5(a));
    Motif_Excitatory_RB_A6(i) = mean(Neuron_RBA6(a));
    Motif_Excitatory_RB_A6_err(i) = var(Neuron_RBA6(a));
    
    Motif_Inhibitory_RB(i) = mean(Neuron_RBA1(b));
    Motif_Inhibitory_RB_err(i) = var(Neuron_RBA1(b));
end

for i = 1:NumMotifs
    DistanceMotifStim(i) = mean(nonzeros(sqrt((MotifX(i)-ElectrodeY(ElectrodeNo)).^2 + (MotifY(i) - ElectrodeX(ElectrodeNo)).^2)));
end
% 
% figure; plot(DistanceMotifStim,Motif_Excitatory_RB_A1,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A2,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A3,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A4,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A5,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A6,'-');
% title('Motif Excitatory RB'); xlabel('Distance Motif To Stimulus'); ylabel('Applied Current AU'); legend('.01','.03','.1','0.125','0.15','.2');


% figure; plot(DistanceMotifStim,Motif_Excitatory_RB_A1-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A2-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A3-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A4-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A5-Motif_Inhibitory_RB,'-');
% hold on; plot(DistanceMotifStim,Motif_Excitatory_RB_A6-Motif_Inhibitory_RB,'-');
% title('Excitatory-Inhibitory RB Difference'); xlabel('Distance Motif To Stimulus'); ylabel('Applied Current AU'); legend('.01','.03','.1','0.125','0.15','.2');


% Options for plotting
options = struct;
options.handle     = figure;
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

Neuron_Excited_Per_Step1 = zeros(numrepeats,length(I0)); % Calculates the number of neurons excited at every step of current
Neuron_Excited_Per_Step2 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step3 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step4 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step5 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step6 = Neuron_Excited_Per_Step1;

for i = 1:length(I0)
    Neuron_Excited_Per_Step1(:,i) = sum(Neuron_RBA1(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step2(:,i) = sum(Neuron_RBA2(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step3(:,i) = sum(Neuron_RBA3(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step4(:,i) = sum(Neuron_RBA4(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step5(:,i) = sum(Neuron_RBA5(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
    Neuron_Excited_Per_Step6(:,i) = sum(Neuron_RBA6(:,Neuron_Excitatory)     <=I0(i),2)*100/length(Neuron_Excitatory);
end

options.color_area = [255 0 0]./255;
options.color_line = [255 0 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step1,options);
hold on
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step2,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step3,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step4,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step5,options);
hold on
options.color_area = [0 0 128]./255;
options.color_line = [0 0 128]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step6,options);

xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Neurons');
legend('delete','0.01','delete','0.03','delete','0.1','delete','0.125','delete','0.15','delete','0.2');


%% Neuron RB Varying with adding more Neurons
%load('B0.mat'); 
load('B1.mat'); load('B2.mat'); load('B3.mat'); load('B4.mat'); load('B5.mat');
Neuron_Excitatory1 = find(repmat([0 1 1 1 1 1],1,NumMotifs) == 1);
Neuron_Excitatory2 = find(repmat([0 0 1 1 1 1 1],1,NumMotifs) == 1);
Neuron_Excitatory3 = find(repmat([0 0 0 1 1 1 1 1],1,NumMotifs) == 1);
Neuron_Excitatory4 = find(repmat([0 0 0 0 1 1 1 1 1],1,NumMotifs) == 1);
Neuron_Excitatory5 = find(repmat([0 0 0 0 0 1 1 1 1 1],1,NumMotifs) == 1);

% Options for plotting
options = struct;
options.handle     = figure;
options.color_area = [0 128 0]./255; % Green : All Neurons
options.color_line = [0 128 0]./255;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

Neuron_Excited_Per_Step1 = zeros(numrepeats,length(I0)); % Calculates the number of neurons excited at every step of current
Neuron_Excited_Per_Step2 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step3 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step4 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step5 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step6 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step7 = Neuron_Excited_Per_Step1; Neuron_Excited_Per_Step8 = Neuron_Excited_Per_Step1;

for i = 1:length(I0)
    %Neuron_Excited_Per_Step1(:,i) = sum(Neuron_RBB0                           <=I0(i),2)*100/length(Neuron_RBB0);
    Neuron_Excited_Per_Step2(:,i) = sum(Neuron_RBB1(:,Neuron_Excitatory1)     <=I0(i),2)*100/length(Neuron_Excitatory1);
    Neuron_Excited_Per_Step3(:,i) = sum(Neuron_RBB2(:,Neuron_Excitatory2)     <=I0(i),2)*100/length(Neuron_Excitatory2);
    Neuron_Excited_Per_Step4(:,i) = sum(Neuron_RBB3(:,Neuron_Excitatory3)     <=I0(i),2)*100/length(Neuron_Excitatory3);
    Neuron_Excited_Per_Step5(:,i) = sum(Neuron_RBB4(:,Neuron_Excitatory4)     <=I0(i),2)*100/length(Neuron_Excitatory4);
    Neuron_Excited_Per_Step6(:,i) = sum(Neuron_RBB5(:,Neuron_Excitatory5)     <=I0(i),2)*100/length(Neuron_Excitatory5);
end

options.color_area = [255 0 0]./255;
options.color_line = [255 0 0]./255;
%plot_areaerrorbar(Neuron_Excited_Per_Step1,options);
hold on
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step2,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step3,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step4,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step5,options);
hold on
options.color_area = [0 0 128]./255;
options.color_line = [0 0 128]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step6,options);
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Excitatory Neurons at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Excitatory Neurons');
legend('delete','1 Inhibitory','delete','2 Inhibitory','delete','3 Inhibitory','delete','4 Inhibitory','delete','5 Inhibitory');

% Subtract 5-4-3-2-1 for 4 total plots
options.handle     = figure;
options.color_area = [0 255 0]./255;
options.color_line = [0 255 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step2-Neuron_Excited_Per_Step3,options);
hold on
options.color_area = [0 0 255]./255;
options.color_line = [0 0 255]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step3-Neuron_Excited_Per_Step4,options);
hold on
options.color_area = [128 0 0]./255;
options.color_line = [128 0 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step4-Neuron_Excited_Per_Step5,options);
hold on
options.color_area = [0 128 0]./255;
options.color_line = [0 128 0]./255;
plot_areaerrorbar(Neuron_Excited_Per_Step5-Neuron_Excited_Per_Step6,options);
xt = get(gca, 'XTick'); set(gca, 'XTick',xt, 'XTickLabel',xt*h); title('Activated Excitatory Difference at Applied Currents'); xlabel('Applied Current AU'); ylabel('Percentage Activated Excitatory Neurons');
legend('delete','1-2 Inhibitory','delete','2-3 Inhibitory','delete','3-4 Inhibitory','delete','4-5 Inhibitory','delete','5-6 Inhibitory');
ylim([0 100]);

%% Functions

function plot_areaerrorbar(data, options)
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