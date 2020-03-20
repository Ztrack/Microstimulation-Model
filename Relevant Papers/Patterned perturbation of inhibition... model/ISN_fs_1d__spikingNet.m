
% - simulating specific ISNs
% in spiking networks of neurons tuned to a 1d feature (orientation)
% patterned inhibitory perturbations delivered based on feature similarity

%% - network params

dt = 1;
t_end = 3e3;
T = 0:dt:t_end;
tau = 20;
vth = 20;

t0 = 1000; tf = 2000;
t_pert = logical((T>t0).*(T<tf));

t_trans = 300;
t_pert_rec = logical((T>(t0+t_trans)).*(T<tf));
t_base_rec = logical((T>(t_trans)).*(T<t0));

NE = 500;
NI = 500;
N = NE+NI;

J0 = 2;
g = -2;
JEE = J0;
JEI = J0;
JIE = g*J0;
JII = g*J0;

mEE = 1;
mEI = 1;
mIE = 1;
mII = 1;

make_weightMatrix = 1;
simulate_network = 1;

%% - weight matrix
if make_weightMatrix

    po_exc = linspace(0,pi,NE);
    po_exc = po_exc(randperm(NE));
    po_inh = linspace(0,pi,NI);
    po_inh = po_inh(randperm(NI));

    po_all = [po_exc, po_inh];

    wEE = zeros(NE,NE);
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

    w = w .* (1+.5*(rand(N)-.5)*2);

    w(eye(N)==1) = 0;

end

%% - activity dynamics

if simulate_network

    b_rate = 2000;

    N_itr = 10;
    s_all = [];
    for itr = 1:N_itr
        itr

        I = poissrnd(b_rate*dt/1000.,N,length(T));

        % - perturbing the specific mode
        p_rate = b_rate*(1-.1*((sin(2*po_inh)')+1));

        for i = 1:NI
            I(NE+i,t_pert) = poissrnd(p_rate(i)*dt/1000,1,sum(t_pert==1));
        end

        [v, s] = LIF_net(I, w, dt, tau, vth);

        s_all = cat(3, s_all, s);
    end
    
end

%% 
function [v,s] = LIF_net(I, W, dt, tau, vth)

    A = -1./tau; %coefficient matrix of the subthreshold dynamics
    B = exp(A*dt);

    s = zeros(size(I));      
    v = zeros(size(I));
    
    for i = 2:size(I,2)
        v(:,i) = B*v(:,i-1) + I(:,i) + (s(:,i-1)'*W)';
        s(:,i) = (v(:,i) > vth);
        v(:,i) = v(:,i) .* (v(:,i) <= vth);
    end
    
end