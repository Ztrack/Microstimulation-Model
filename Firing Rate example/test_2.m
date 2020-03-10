clearvars; close all;
num_stim = 5; % Number of stimulations in present in matrix area
sx= 100;
sy= 100;

stim_loc =  zeros(sx, sy); % Initialize stimulation location matrix size nxm
re = randperm(numel(stim_loc),num_stim); % Selects random element n times
stim_loc(re) = 1;

[axon_x,axon_y] = find(stim_loc);
axon_loc =  zeros(sx, sy);
for i = 1:length(axon_x)
    x = [axon_x(i) randi([1,100],1)];                   % x coordinates
    y = [axon_y(i) randi([1,100],1)];                   % y coordinates
    Points = max(abs(diff(x)), abs(diff(y)))+1;    % Number of points in line
    rIndex = round(linspace(y(1), y(2), Points));  % Row indices
    cIndex = round(linspace(x(1), x(2), Points));  % Column indices
    lindex = sub2ind(size(axon_loc), rIndex, cIndex);     % Linear indices
    axon_loc(lindex) = 1;  % Set the line points to white
end


figure; imagesc(axon_loc); title('Excited Axons'); xlabel('X Position'); ylabel('Y Position')