
data = load('solrbhomo.mat');

load('InitialConditionsFull.mat');

clearvars -except neuron data;

c_mat = data.solrbhomo.c;
c_units = data.solrbhomo.units_c;
o_mat = data.solrbhomo.o;
o_units = data.solrbhomo.units_o;


clear data;

neu_type = neuron.adapt.type;
rf = neuron.adapt.ratefunction;

%% plot (current required for 50% activation)

current_steps = [0.04:0.04:2];

excit_vals = zeros(1,10);
for n=1:10
    a = mean(c_mat(:,neuron.excitatory,n),2);
    for m=1:length(a)
        if a(m) > 50
            excit_vals(n) = m;
            break;
        end
    end
end

inhib_vals = zeros(1,10);
for n=1:10
    a = mean(c_mat(:,neuron.inhibitory,n),2);
    for m=1:length(a)
        if a(m) > 50
            inhib_vals(n) = m;
            break;
        end
    end
end

figure; %FIGURE 1
% subplot(2,2,1);
e = bar([100:100:1000],current_steps(excit_vals),'blue'); hold on;
i = bar([100:100:1000],current_steps(inhib_vals),'red');
legend([e i],{'Excitatory','Inhibitory'});
ylabel('Current step (mA) x10^5');
xlabel('Bin # (100ms)');
title('Current required for 50% activation');
hold off;


%% other plots (current / excitatory)

subplot(2,2,1);
l1 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,1),2),'LineWidth',2,'Color',[1 0 0]); hold on; 
l2 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,2),2),'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,3),2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,4),2),'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,5),2),'LineWidth',2,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,6),2),'LineWidth',2,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,7),2),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,8),2),'LineWidth',2,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,9),2),'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,neuron.excitatory,10),2),'LineWidth',2,'Color',[1 0 1]);

legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');

title('Excitatory Neuron Activation Through the Center Electrode', 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;

%% rate function plots

for n=1:3
    rf(n,:) = (rf(n,:) .* 0.9) + 3; %shift up to add baseline
end

figure; %FIGURE 2

subplot(1,3,1);
ra = plot([1:1000],rf(1,:),'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
title('Rapidly Adapting Function (562 units)', 'FontSize', 12);
ylabel('Spikes per Second');
xlabel('Time (ms)');

subplot(1,3,2);
sa = plot([1:1000],rf(2,:),'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
title('Slowly Adapting Function (38 units)', 'FontSize', 12);
ylabel('Spikes per Second');
xlabel('Time (ms)');

subplot(1,3,3);
mi = plot([1:1000],rf(3,:),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
title('Rapidly Adapting Function (400 units)', 'FontSize', 12);
ylabel('Spikes per Second');
xlabel('Time (ms)');


% plot([1:100],rateFunction(1:100),'LineWidth',2,'Color',[1 0 0]); hold on;
% plot([101:200],rateFunction(101:200),'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
% plot([201:300],rateFunction(201:300),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
% plot([301:400],rateFunction(301:400),'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
% plot([401:500],rateFunction(401:500),'LineWidth',2,'Color',[1 1 0]);
% plot([501:600],rateFunction(501:600),'LineWidth',2,'Color',[0 1 0]);
% plot([601:700],rateFunction(601:700),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
% plot([701:800],rateFunction(701:800),'LineWidth',2,'Color',[0 0.4470 0.7410]);
% plot([801:900],rateFunction(801:900),'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
% plot([901:1000],rateFunction(901:1000),'LineWidth',2,'Color',[1 0 1]);
% title('Non-Homogeneous Rate Function','FontSize', 14);
% xlabel('Time (ms)');
% ylabel('Spikes / sec'); yticks([0:5:35]);


%% current / inhibitory

subplot(2,2,2);
l1 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,1),2),'LineWidth',2,'Color',[1 0 0]); hold on; 
l2 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,2),2),'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,3),2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,4),2),'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,5),2),'LineWidth',2,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,6),2),'LineWidth',2,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,7),2),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,8),2),'LineWidth',2,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,9),2),'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,neuron.inhibitory,10),2),'LineWidth',2,'Color',[1 0 1]);

legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');

title('Inhibitory Neuron Activation Through the Center Electrode', 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;


%% RATE AND SUBTYPE PLOTS
% 20200806 - with standard error maggin shading

options=struct();
options.handle     = figure; %FIGURE 3
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
options.x_axis     = [0.04:0.04:2];

colors = [0.6350 0.0780 0.1840; 1 0 0; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    1 1 0; 0 1 0; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.4940 0.1840 0.5560;
    1 0 1];

% l1 = plot([0.04:0.04:2],mean(c_mat(:,:,1),2),'LineWidth',2,'Color',[0.6350 0.0780 0.1840]); hold on;
% l2 = plot([0.04:0.04:2],mean(c_mat(:,:,2),2),'LineWidth',2,'Color',[1 0 0]);
% l3 = plot([0.04:0.04:2],mean(c_mat(:,:,3),2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
% l4 = plot([0.04:0.04:2],mean(c_mat(:,:,4),2),'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
% l5 = plot([0.04:0.04:2],mean(c_mat(:,:,5),2),'LineWidth',2,'Color',[1 1 0]);
% l6 = plot([0.04:0.04:2],mean(c_mat(:,:,6),2),'LineWidth',2,'Color',[0 1 0]);
% l7 = plot([0.04:0.04:2],mean(c_mat(:,:,7),2),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
% l8 = plot([0.04:0.04:2],mean(c_mat(:,:,8),2),'LineWidth',2,'Color',[0 0.4470 0.7410]);
% l9 = plot([0.04:0.04:2],mean(c_mat(:,:,9),2),'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
% l0 = plot([0.04:0.04:2],mean(c_mat(:,:,10),2),'LineWidth',2,'Color',[1 0 1]);

options.color_area = [0.6350 0.0780 0.1840]; options.color_line = [0.6350 0.0780 0.1840];
l1 = plot_areaerrorbar(c_mat(:,:,1)',options); hold on;
options.color_area = [1 0 0]; options.color_line = [1 0 0];
l2 = plot_areaerrorbar(c_mat(:,:,2)',options);
options.color_area = [0.8500 0.3250 0.0980]; options.color_line = [0.8500 0.3250 0.0980];
l3 = plot_areaerrorbar(c_mat(:,:,3)',options);
options.color_area = [0.9290 0.6940 0.1250]; options.color_line = [0.9290 0.6940 0.1250];
l4 = plot_areaerrorbar(c_mat(:,:,4)',options);
options.color_area = [1 1 0]; options.color_line = [1 1 0];
l5 = plot_areaerrorbar(c_mat(:,:,5)',options);
options.color_area = [0 1 0]; options.color_line = [0 1 0];
l6 = plot_areaerrorbar(c_mat(:,:,6)',options);
options.color_area = [0.4660 0.6740 0.1880]; options.color_line = [0.4660 0.6740 0.1880];
l7 = plot_areaerrorbar(c_mat(:,:,7)',options);
options.color_area = [0 0.4470 0.7410]; options.color_line = [0 0.4470 0.7410];
l8 = plot_areaerrorbar(c_mat(:,:,8)',options);
options.color_area = [0.4940 0.1840 0.5560]; options.color_line = [0.4940 0.1840 0.5560];
l9 = plot_areaerrorbar(c_mat(:,:,9)',options);
options.color_area = [1 0 1]; options.color_line = [1 0 1];
l0 = plot_areaerrorbar(c_mat(:,:,10)',options);


%legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
%    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
%    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');

title(' (All) Neuron Activation Through the Center Electrode', 'FontSize', 12);
ylabel('Percent Activation');
xlabel(strcat("Current (nA / mm^2",")    x 10^5")); xticks([0 0.5 1 1.5 2.0]);
hold off;

%%  --------------------------------------------------------------------------------------------
figure; %FIGURE 4

% neuron.type==1 is Inhibitory
% neuron.type==2 is Excitatory
%
% neu_type == 1 is RA
% neu_type == 2 is SA
% neu_type == 3 is Mixed


subplot(2,3,1);
cond = neuron.type==1 & neu_type==1;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Inhibitory RA Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subplot(2,3,2);
cond = neuron.type==1 & neu_type==2;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Inhibitory SA Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subplot(2,3,3);
cond = neuron.type==1 & neu_type==3;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Inhibitory Mixed Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subplot(2,3,4);
cond = neuron.type==2 & neu_type==1;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Excitatory RA Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subplot(2,3,5);
cond = neuron.type==2 & neu_type==2;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Excitatory SA Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subplot(2,3,6);
cond = neuron.type==2 & neu_type==3;
l1 = plot([0.04:0.04:2],mean(c_mat(:,cond,1),2),'LineWidth',1,'Color',[0.6350 0.0780 0.1840]); hold on;
l2 = plot([0.04:0.04:2],mean(c_mat(:,cond,2),2),'LineWidth',1,'Color',[1 0 0]);
l3 = plot([0.04:0.04:2],mean(c_mat(:,cond,3),2),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
l4 = plot([0.04:0.04:2],mean(c_mat(:,cond,4),2),'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
l5 = plot([0.04:0.04:2],mean(c_mat(:,cond,5),2),'LineWidth',1,'Color',[1 1 0]);
l6 = plot([0.04:0.04:2],mean(c_mat(:,cond,6),2),'LineWidth',1,'Color',[0 1 0]);
l7 = plot([0.04:0.04:2],mean(c_mat(:,cond,7),2),'LineWidth',1,'Color',[0.4660 0.6740 0.1880]);
l8 = plot([0.04:0.04:2],mean(c_mat(:,cond,8),2),'LineWidth',1,'Color',[0 0.4470 0.7410]);
l9 = plot([0.04:0.04:2],mean(c_mat(:,cond,9),2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
l0 = plot([0.04:0.04:2],mean(c_mat(:,cond,10),2),'LineWidth',1,'Color',[1 0 1]);
legend([l1 l2 l3 l4 l5 l6 l7 l8 l9 l0],{'(1-100ms)','(101-200ms)',...
    '(201-300ms)','(301-400ms)','(401-500ms)','(501-600ms)',...
    '(601-700ms)','(701-800ms)','(801-900ms)','(901-1000ms)'},'Location','southeast');
title(strcat("Excitatory Mixed Neurons ",num2str(sum(cond))," units"), 'FontSize', 12);
ylabel('Percent Activation');
xlabel('Current (mA)    x 10^5'); xticks([0 0.5 1 1.5 2.0]);
hold off;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 






