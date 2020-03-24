
% - simulating specific ISNs
% in rate-based networks of neurons tuned to a 1d feature (orientation)
% patterned inhibitory perturbations delivered based on feature similarity

%% - network params

dt = 1; % 1ms
t_end = 1e3; % 1000ms
T = 0:dt:t_end; % Time vector
tau = 10; % 

t0 = 300; tf = 700;
t_pert = logical((T>t0).*(T<tf)); % Creates array of 1's between time 300 and 700

t_trans = 100;
t_pert_rec = logical((T>(t0+t_trans)).*(T<tf));  % ARray of 1's between 400:700
t_base_rec = logical((T>(t_trans)).*(T<t0)); % ARray of 1's between 100:300

NE = 400; % Neuron excitatory
NI = 400; % Neuron inhibitory
N = NE+NI; % Total neurons

J0 = .1/2;
JEE = J0; % Patterned perturbations?
JEI = J0; 
JIE = -1.5*J0; %  Randomized patterns over neurons?
JII = -1.5*J0;

mEE = 1; % Parameters of connectivity ???
mEI = 1;
mIE = 1;
mII = 1;

po_exc = linspace(0,pi,NE); % Pi spaced apart from 0 to # Neurons? wtf?
po_inh = linspace(0,pi,NI);

simulate_network = 1;

%% - weight matrix

wEE = zeros(NE,NE); % Tuning curve of the neuron with respect to other neurons from 0:1
wEI = zeros(NE,NI);
for i = 1:NE
    wEE(i,:) = (1 + mEE * cos(2*(po_exc(i) - po_exc))) * JEE;
    wEI(i,:) = (1 + mEI * cos(2*(po_exc(i) - po_inh))) * JEI;
end
wIE = zeros(NI,NE);
wII = zeros(NI,NI);
for i = 1:NI
    wIE(i,:) = (1 + mIE * cos(2*(po_inh(i) - po_exc))) * JIE;
    wII(i,:) = (1 + mII * cos(2*(po_inh(i) - po_inh))) * JII;
end

w = [wEE, wEI
     wIE, wII];
 
w = w + J0/2*(rand(N,N)-.5);

w(1:NE,:) = rectify(w(1:NE,:));
w(NE+1:end,:) = -rectify(-w(NE+1:end,:));

w(eye(N)==1) = 0;

figure; plot(wEE(1,:)); % Visualize curve

%% - activity dynamics

if simulate_network

    % - perturbing the specific mode
    I_specPert = 1 + .25*rand(N,length(T))/5; % N by time long % change in the input to inhibitory (inh.) neurons during pert.
    dI_specPert = -.1*(sin(2*po_inh)')-.1; 
    I_specPert(NE+1:end,t_pert) = I_specPert(NE+1:end,t_pert) + dI_specPert; % input to pert. inh. neurons

    r_specPert = simulate_dynamics(I_specPert, N, T, w, dt, tau); % output activity of pert. inh. neurons

    % - perturbing the uniform mode0
    I_genPert = 1 + .25*rand(N,length(T))/5;
    dI_genPert = dI_specPert(randperm(NI)); % randomizing the specific pert.
    I_genPert(NE+1:end,t_pert) = I_genPert(NE+1:end,t_pert) + dI_genPert;

    r_genPert = simulate_dynamics(I_genPert, N, T, w, dt, tau);

end

%% - functions
function [r] = simulate_dynamics(I, N, T, w, dt, tau)
    r = zeros(N, length(T));
    for i = 1:length(T)-1
        inp = w' * r(:,i) + I(:,i);
        rd = -r(:,i) + rectify(inp);
        r(:,i+1) = r(:,i) + dt/tau * rd;
    end
end

function [zr] = rectify(z)
    zr = z .* (z>0);
end
