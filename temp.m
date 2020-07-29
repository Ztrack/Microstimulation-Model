%% Plot temporary

% Options for plotting
options = struct;
options.alpha      = 0.5;
options.line_width = 2;
options.error = 'c95';
options.legendswitch = 0; % Legend 0 = off, 1 = on
options.legend = [];

options.handle     = figure;
options.color_area = [0 128 0]./255; % Green
plot_areaerrorbar(xsol1curve',options); 
title('Stimulation Performance');
xlabel('Iteration Number');
ylabel('Non-Motion / Motion-Tuned Neuron Ratio');
legend('Microstimulation');
ylim([0 6]);
hold off

options.handle     = figure;
plot_areaerrorbar(xsol1curve',options);
hold on;
options.color_area = [255 0 0]./255; % Red
plot_areaerrorbar(xsol4curve',options);
ylim([0 6]);
hold off
title('Stimulation Performance');
xlabel('Iteration Number');
ylabel('Non-Motion / Motion-Tuned Neuron Ratio');
legend('Microstimulation','Microstimulation + Optogenetics');






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