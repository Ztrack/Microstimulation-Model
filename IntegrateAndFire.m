function out_FR = IntegrateAndFire(neuron,params,lambdahat)

% Spike Model
I = zeros(params.numneurons,1); % Initialize synaptic input from other neurons. 0 at time start.
S = zeros(1,params.numneurons); % Initialize a binary spiking vector of the neurons spiking in previous time step
out_FR = zeros(params.numneurons,1); % Initialize output rate
nbins = 1000; % Time vector length
%raster = zeros(1000,length(S));

% Neuron type data taken from https://www.neuroelectro.org/
tP = zeros(params.numneurons,1); tP(neuron.type == 1) = 9.57; tP(neuron.type == 2) = 19.73; % Membrane Time constant (ms)- Time constant for the membrane to repolarize after a small current injection of fixed amplitude and duration
Pr = zeros(params.numneurons,1); Pr(neuron.type == 1) = 40; Pr(neuron.type == 2) = 20; % Background spiking probability (Hz) - AP discharge rate in the absence of current injection or a stimulus
Pr = (Pr + lambdahat)/nbins;
Ps = Pr; % Initialize probability of spiking

for t = 1:1:nbins
    
    % Each time step (dt) of 1 ms starts with each neuron, i, updating Ii with received input from any of the set of connected neurons, together with an exponential synaptic decay, which can either be excitatory or inhibitory.

     % Update Synaptic Weights
    I = I + sum(neuron.weight.matrix.*(S./nbins),2); % if previous bin spikes, S=1 applying a synaptic effect from that neuron.
    
    % Determine Spiking
    S = Ps + I > rand(params.numneurons,1); % We determine whether the neuron spikes with the probability PS and update the spiking vector for the next time step
    
    % Update params
    dIdt = (0 - I)./tP; % Synaptic decay with tI time factor
    I = I + dIdt; % Apply spiking decay
    I(S == 1) = 0; % Reset synaptic strength if neuron fired
    Ps(S == 1) = Pr(S == 1); %If a neuron spikes, the probability is reset to the reset value
    
    % Output
    %raster(t,:) = S; % For debugging
    out_FR = out_FR + S;
end

%% Optional Plotting - Debugging

% figure;
% plot(1:1000,sum(raster,2));
% title('Population Activity');
% ylabel('Number of Spikes');
% xlabel('Time (ms)');
% 
% figure;
% hist(sum(raster,1),25);
% title('Neuron Activity Across Population');
% ylabel('Number of Neurons');
% xlabel('Firing Rate (Hz)');
% 
% figure;
% for i = 1:params.numneurons
%     ylim([1 200]);
%     xlim([0 1000]);
%     plot(1:1000,raster(:,i).*i,'.');
%     hold on;
% end
% title('Raster Data');
% ylabel('Neuron Number');
% xlabel('Time (ms)');
% 
% figure;
% x = 1:2;
% sumspikes = sum(raster,1);
% y = [mean(sum(raster(:,neuron.type == 1),1)) mean(sum(raster(:,neuron.type == 2),1))];
% ystd = [std(sum(raster(:,neuron.type == 1),1)) std(sum(raster(:,neuron.type == 2),1))];
% y95 = 1.96*[ystd(1)/sqrt(sum(neuron.type==1)) ystd(2)/sqrt(sum(neuron.type==2))];
% bar(x,y); hold on;
% er = errorbar(x,y,ystd,ystd);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% title('Population Firing Rate');
% ylabel('firing rate (Hz)');
% xlabel('Inhibitory                Excitatory');
% ylim([0 max(y)+10]);

end