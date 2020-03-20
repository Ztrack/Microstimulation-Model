
% visualizes the output of patterned perturbations of rate-based networks
% obtained from, e.g., ISN_fs_1d.m

%% - figure params

fig_pos = [100,100,500,200];
fs = 15;
tkl = .025;
cl_exc = 'm';
cl_inh = [0,0.5,.9];

print_res = '-r300';

%% - visualize the dynamics

% - parameters to visualize:
% change in the input to inhibitory (inh.) neurons during pert.
dI_pert = dI_specPert;
% input to pert. inh. neurons
I_pert = I_specPert;
% output activity of pert. inh. neurons
r_pert = r_specPert;

plot_generalActivity(r_pert, T, NE, fig_pos, cl_exc, cl_inh, fs, tkl);
print(['GeneralActivity_SpecPert__' num2str(mEE)], '-dpng', print_res)

plot_specActivity(r_pert, t_pert_rec, t_base_rec, I_pert, po_exc, po_inh, NE, fig_pos, cl_exc, cl_inh, fs, tkl);
print(['SpecActivity_SpecPert__' num2str(mEE)], '-dpng', print_res)

plot_regress_line = 1;
plot_pertChange(dI_pert, r_pert, t_pert_rec, t_base_rec, NE, fig_pos, fs, tkl, cl_inh, plot_regress_line);
ylim([-.1,.4])
print(['PertChange_SpecPert__' num2str(mEE)], '-dpng', print_res)


%% - functions
function [z] = plot_generalActivity(r, T, NE, fig_pos, cl_exc, cl_inh, fs, tkl)
z = [];

figure('Position', fig_pos);

subplot(121); title('Exc', 'color', cl_exc); hold on
plot(T,r(1:NE,:), 'color', cl_exc, 'LineWidth',.25)
plot(T,nanmean(r(1:NE,:),1), 'k', 'LineWidth',2)
xlabel('Time (a.u.)')
ylabel('Activity (a.u.)')
ylim([0,.6])
yticks([0,.5])

fig_format(fs,tkl)

subplot(122); title('Inh', 'color', cl_inh); hold on
plot(T,r(NE+1:end,:), 'color', cl_inh, 'LineWidth',.25)
plot(T,nanmean(r(NE+1:end,:),1), 'k', 'LineWidth',2)
xlabel('Time (a.u.)')
ylim([0,.6])
yticks([0,.5])

fig_format(fs,tkl)

end

function [z] = plot_pertChange(dI_pert, r_pert, t_pert_rec, t_base_rec, NE, fig_pos, fs, tkl, cl_inh, plot_regress_line)
    figure('Position', [100,100,fig_pos(3)/2,fig_pos(4)]);
    subplot(111); title('Inh', 'color', cl_inh); 
    hold on;
    
    x = dI_pert;
    x0 = min(x(:))-.01; xf = max(x(:))+.01;
    xx = x0:.01:xf;
    
    y = nanmean(r_pert(NE+1:end,t_pert_rec),2) - nanmean(r_pert(NE+1:end,t_base_rec),2);
    y0 = min(y(:))-.1; yf = max(y(:))+.1;
    
    plot(x, y, 'o', 'color', cl_inh)
    
    if plot_regress_line
        X = [ones(length(x),1) x];
        r = X\y;

        plot(xx, r(1) + r(2)*xx, 'r-', 'LineWidth',2);
        text(0.1,1,num2str(round(r(2),1)), 'FontSize',fs,'Units', 'normalized', 'Color','r')
    end
    
    xlabel('Input perturbation')
    ylabel('Response change')
    
    xlim([x0,xf])
    
    fig_format(fs,tkl)

end

function [z] = plot_specActivity(r, t_pert_rec, t_base_rec, I, po_exc, po_inh, NE, fig_pos, cl_exc, cl_inh, fs, tkl)
z = [];

figure('Position', fig_pos);

subplot(121); title('Exc', 'color', cl_exc); 
hold on
plot(po_exc*180/pi, nanmean(r(1:NE,t_base_rec),2), 'o', 'color', 'k')
plot(po_exc*180/pi, nanmean(r(1:NE,t_pert_rec),2), 'color', cl_exc)
xlabel('Exc. pref. orient. (deg)')
ylabel('Avg. activity (a.u.)')
xticks([0,45,90,135,180])
xticklabels({'0','','90','','180'})
ylim([0,.7])
xlim([0,180])
legend('Baseline', 'Perturb')
legend boxoff

fig_format(fs,tkl)

subplot(122); title('Inh', 'color', cl_inh); 
hold on
plot(po_inh*180/pi, nanmean(r(NE+1:end,t_base_rec),2), 'o', 'color', 'k')
plot(po_inh*180/pi, nanmean(r(NE+1:end,t_pert_rec),2), 'color', cl_inh)
xlabel('Inh. pref. orient. (deg)')
xticks([0,45,90,135,180])
xticklabels({'0','','90','','180'})
ylim([0,.7])
xlim([0,180])

yyaxis right; hold on
plot(po_inh*180/pi, nanmean(I(NE+1:end,t_pert_rec),2)-nanmean(I(NE+1:end,t_base_rec),2), ...
    'LineWidth',1)
ylim([-.5,0])
yticks([-.2,-.1,0])
yticklabels({-.2,'',0})
ylabel('Perturbation')

fig_format(fs,tkl)
end

function [] = fig_format(fs, tkl)
    set(gca, 'LineWidth', 1, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'TickLength', [tkl tkl])
end