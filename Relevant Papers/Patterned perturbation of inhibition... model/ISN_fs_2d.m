
% - simulating specific ISNs
% in rate-based networks of neurons with 2d visual RFs  
% patterned inhibitory perturbations delivered based on feature similarity (1d or 2d)

%% - network params

dt = 1;
t_end = 1e3;
T = 0:dt:t_end;
tau = 10;

t0 = 300; tf = 700;
t_pert = logical((T>t0).*(T<tf));

t_trans = 100;
t_pert_rec = logical((T>(t0+t_trans)).*(T<tf));
t_base_rec = logical((T>(t_trans)).*(T<t0));

NE = 400;
NI = 400;
N = NE+NI;

J0 = .1/2;
JEE = J0;
JEI = J0;
JIE = -1.5*J0;
JII = -1.5*J0;

mEE = .5;
mEI = .5;
mIE = .5;
mII = .5;

exp_pwr = 2;

%

make_RFs = 1;
make_weigthMatrix = 1;
simulate_network = 1;
    
%% - receptive fields
if make_RFs
    
    ppd = 4; % pixel per degree
    vf_size = 50; % in degrees (visual field)
    sz = vf_size*ppd; % size in pixel

    locs = (rand(N, 2)-.5)*(sz/20) ;

    po_exc = linspace(0,pi,NE);
    po_inh = linspace(0,pi,NI);
    po_all = [po_exc, po_inh];

    sigmas = 2.5*ones([1, N]);
    gammas = .5*ones([1, N]);

    psis = rand(1,N) * pi;       
    sfs = gamrnd(2, 1, [1, N]) * .04; 
    sfs(sfs < .01) = .01;
    sfs(sfs > 1) = 1;

    RFs = zeros(N, sz, sz);
    for i = 1:N
        rf = GaborRFs(sigmas(i), gammas(i), psis(i), 1/sfs(i), po_all(i), locs(i,:), ppd, vf_size);
        RFs(i,:,:) = rf(1:sz,1:sz);
    end

end

%% - weight matrix
if make_weigthMatrix
    
    cc_rfs = corr(reshape(RFs, N, [])');
    cc_rfs(eye(N)==1)=nan;

    wEE = JEE* (.1 + mEE* exp(exp_pwr* cc_rfs(1:NE,1:NE)));
    wEI = JEI* (.1 + mEI* exp(exp_pwr* cc_rfs(1:NE,NE+1:end)));
    wIE = JIE* (.1 + mIE* exp(exp_pwr* cc_rfs(NE+1:end,1:NE)));
    wII = JII* (.1 + mII* exp(exp_pwr* cc_rfs(NE+1:end,NE+1:end)));

    w = [wEE, wEI
         wIE, wII];

    w = w + J0/2*(rand(N,N)-.5);

    w(1:NE,:) = rectify(w(1:NE,:));
    w(NE+1:end,:) = -rectify(-w(NE+1:end,:));

    w(eye(N)==1) = 0;

    cc_rfs(eye(N)==1)=0;

end


%% - activity dynamics
if simulate_network

    % - perturbing the specific mode - along 1d orientation
    I_specPert = 1 + .25*rand(N,length(T))/5;
    dI_pert_1d = -.2*sin(2*po_inh)' -.2;

    dI_pert_1d(dI_pert_1d>0)=0;
    I_specPert(NE+1:end,t_pert) = I_specPert(NE+1:end,t_pert) + dI_pert_1d;

    r_specPert_1d = simulate_dynamics(I_specPert, N, T, w, dt, tau);

    
    % - perturbing the specific mode - along 2d RFs
    nid = NE+2;
    inh_dim = cc_rfs(nid,NE+1:end);

    I_specPert = 1 + .25*rand(N,length(T))/5;
    dI_pert = -.2*(exp(1*inh_dim'));
    I_specPert(NE+1:end,t_pert) = I_specPert(NE+1:end,t_pert) + dI_pert;
    I_specPert(I_specPert<0) = 0;

    r_specPert = simulate_dynamics(I_specPert, N, T, w, dt, tau);

    
    % - perturbing the specific mode - along 2d RFs (randomized control)
    inh_dim = inh_dim(randperm(length(inh_dim)));

    I_specPert_randomized = 1 + .25*rand(N,length(T))/5;
    dI_pert = -.2*(exp(1*inh_dim'));
    I_specPert_randomized(NE+1:end,t_pert) = I_specPert_randomized(NE+1:end,t_pert) + dI_pert;
    I_specPert_randomized(I_specPert<0) = 0;

    r_specPert_randomized = simulate_dynamics(I_specPert_randomized, N, T, w, dt, tau);

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
