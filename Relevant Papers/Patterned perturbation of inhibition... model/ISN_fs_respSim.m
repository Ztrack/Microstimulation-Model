
% - simulating specific ISNs
% in rate-based networks of neurons with 2d visual RFs 
% patterned inhibitory perturbations delivered based on response similarity

%% - network params

dt = 1;
t_end = 1e3;
T = 0:dt:t_end;
tau = 10;

t0 = 500; tf = 1000;
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

n_stim = 200;
t_each = 200;

%

make_RFs = 1;
make_weigthMatrix = 1;
make_stims = 1;
response_correlations = 1;
perturb_network_sample = 1;

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

%%
if make_stims
    
    locs_stim = zeros(n_stim, 2);

    po_stim = linspace(0,pi,n_stim);
    po_stim = po_stim(randperm(n_stim));

    sigmas_stim = 2.5*ones([1, n_stim]);
    gammas_stim = .5*ones([1, n_stim]);
    psis_stim = rand(1,n_stim) * pi;         
    sfs_stim = gamrnd(2, 1, [1, n_stim]) * .04; 

    stims = zeros(n_stim, sz, sz);
    for i = 1:n_stim
        rf = GaborRFs(sigmas_stim(i), gammas_stim(i), psis_stim(i), 1/sfs_stim(i), po_stim(i), locs_stim(i,:), ppd, vf_size);
        stims(i,:,:) = rf(1:sz,1:sz);

        % blank, every other one
        if mod(i,2) == 0
           stims(i,:,:) = nanmean(nanmean(stims(i,:,:)));
        end
    end

    cc_sts = corr(reshape(stims, n_stim, [])', reshape(RFs, N, [])'); 

end

%% - response cc
if response_correlations
        
    T_resp = 0:dt:t_each*n_stim-dt;

    I_b = [];
    for i = 1:n_stim
        I_b0 = 1+2*cc_sts(i,:)' + .025*rand(N,t_each)/2;
        I_b = [I_b, I_b0];
    end

    r_b = simulate_dynamics(I_b, N, T_resp, w, dt, tau);

    re_th = prctile(nanmean(r_b(1:NE,:),2),20);
    resp_ids_exc = find(nanmean(r_b(1:NE,:),2)>re_th);
    ri_th = prctile(nanmean(r_b(NE+1:end,:),2),20);
    resp_ids_inh = NE+find(nanmean(r_b(NE+1:end,:),2)>ri_th);

    cc_resp = corr(r_b');

    ccm_resp_exc = corr(reshape(cc_rfs(resp_ids_exc,resp_ids_exc),1,[])', ...
                        reshape(cc_resp(resp_ids_exc,resp_ids_exc),1,[])');
    ccm_resp_inh = corr(reshape(cc_rfs(resp_ids_inh,resp_ids_inh),1,[])', ...
                        reshape(cc_resp(resp_ids_inh,resp_ids_inh),1,[])');
                        

end

%% - activity dynamics (sample perturb.)
if perturb_network_sample
    [a,b] = sort(nanmean(abs(cc_resp(NE+1:end,NE+1:end))));
    nid = NE+b(end-40);
    
    inh_dim = cc_resp(nid,NE+1:end);
    exc_dim = cc_resp(nid,1:NE);

    I_specPert = 1 + .25*rand(N,length(T))/5;
    dI_pert_resp = -.25*(exp(2*inh_dim'));
    I_specPert(NE+1:end,t_pert) = I_specPert(NE+1:end,t_pert) + dI_pert_resp;
    I_specPert(I_specPert<0) = 0;

    r_specPert_resp = simulate_dynamics(I_specPert, N, T, w, dt, tau);

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
