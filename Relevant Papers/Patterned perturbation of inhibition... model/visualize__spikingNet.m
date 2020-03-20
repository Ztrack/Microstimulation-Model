
% visualizes the output of the patterned perturbations of spiking networks
% obtained from: ISN_fs_1d__spikingNet.m

%% - figure params

fig_pos1 = [100,100,800,400];
fig_pos2 = [100,100,350,250];
fig_pos3 = [100,100,250,250];

fs = 15;
tkl = .025;
cl_exc = 'm';
cl_inh = [0,0.5,.9];

res = '-r300';

%% - average rates

r_exc_base = 0; r_exc_pert = 0;
r_inh_base = 0; r_inh_pert = 0;
inp_pert = 0;
for itr = 1:N_itr
    s = squeeze(s_all(:,:,itr));
    
    r_exc_base = r_exc_base+nansum(s(1:NE,t_base_rec),2)/(nansum(t_base_rec)*dt/1e3);
    r_exc_pert = r_exc_pert+nansum(s(1:NE,t_pert_rec),2)/(nansum(t_pert_rec)*dt/1e3);

    r_inh_base = r_inh_base+nansum(s(NE+1:end,t_base_rec),2)/(nansum(t_base_rec)*dt/1e3);
    r_inh_pert = r_inh_pert+nansum(s(NE+1:end,t_pert_rec),2)/(nansum(t_pert_rec)*dt/1e3);

    inp_pert = inp_pert+nanmean(I(NE+1:end,t_pert_rec),2)-nanmean(I(NE+1:end,t_base_rec),2);
end

%% - sample spiking activity

[si,st] = find(s == 1);

figure('Position', fig_pos1)
hold on
plot(st(si<=NE)*dt/1e3,si(si<=NE),'.', 'color', cl_exc)
plot(st(si>NE)*dt/1e3,si(si>NE),'.', 'color', cl_inh)

plot(st(si>N)*dt/1e3,si(si>N),'.', 'color', 'k')

text(-.15,750, 'Inh', 'Color', cl_inh, 'FontSize',fs)
text(-.15,250, 'Exc', 'Color', cl_exc, 'FontSize',fs)

yticks([1, NE, N])

xlabel('Time (s)')
ylabel('Neuron #')

xlim([0, T(end)/1e3])

set(gca, 'LineWidth', 1, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01])

print(['SpikingActivity.png'], '-dpng', '-r300')

%% - plot specific activity
figure('Position', fig_pos2);

subplot(111); title('Inh'); hold on
plot(po_inh*180/pi, (r_inh_pert-r_inh_base)/N_itr, 'o', 'color', cl_inh)
ylabel({'Avg. rate change', '(spikes/s)'})
xlabel('Pref. orient. (deg)')
xticks([0,45,90,135,180])
xticklabels({'0','','90','','180'})
xlim([0,180])

yyaxis right; hold on
plot(po_inh*180/pi, p_rate - b_rate, '.')
ylabel({'Input perturbation', '(spikes/s)'})

set(gca, 'LineWidth', 1, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'TickLength', [tkl tkl])

print(['SpikingNet_SpecActivity_SpecPert__inh.png'], '-dpng', res)

%% - diff. pert. vs diff. activity change (tuning curve)

figure('Position', fig_pos3); 
subplot(111); title('Inh'); hold on;

x = p_rate - b_rate;

x0 = min(x(:))-.01; xf = max(x(:))+.01;
xx = x0:.01:xf;

y = (r_inh_pert - r_inh_base)/N_itr;

y0 = min(y(:))-.1; yf = max(y(:))+.1;

plot(x, y, 'ko')

X = [ones(length(x),1) x];
r = X\y;

plot(xx, r(1) + r(2)*xx, 'r-', 'LineWidth',2);
%text(0.1,1,num2str(round(r(2),2)), 'FontSize',fs,'Units', 'normalized', 'Color','r')

xlabel({'Input perturbation', '(spikes/s)'})
ylabel({'Avg. rate change' ,'(spikes/s)'})

xlim([x0,xf])

set(gca, 'LineWidth', 1, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'TickLength', [tkl tkl])

print(['SpikingNet_PertChange_SpecPert.png'], '-dpng', res)
