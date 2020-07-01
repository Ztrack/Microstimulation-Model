function output = thresholdfun(SEoutput,units,n,x)
% Function to convert raw output to thresholds
% output = num electrodes x num steps x num neurons x num simulations
% output = 100 x 50 x 1000 x 2000 in typical settings

%% Constructing theshold per electrode
% We must find a proper threshold for each elecrode. This will be based
% on the first x neurons of the motion orientation we are interested in
% new output must be a low and a high value per stimulation strategy (4 total)
% n = 0.5 % Percent activation to count as 'activated'
% x = The current/opto step must activate this many motion tuned neurons to store the step

output = zeros(1,size(SEoutput,1));

for i = 1:size(SEoutput,1) % For every electrode
    numactivated = zeros(1,size(SEoutput,2));
    
    for ii = 1:size(SEoutput,2) % For every step
        activated = squeeze(SEoutput(i,ii,:) > max(max(SEoutput))*n);
        numactivated(ii) = sum(activated) > x;
    end
    
    output(i) = units(find(numactivated == 1,1));
end



%%
% Calculating distribution of 50% threshold for all neurons
% dist1 = zeros(size(SEoutput));
% for i = 1:size(SEoutput,1) % For every electrode
%     
%     for ii = 1:size(SEoutput,2) % For every step
%         activated = squeeze(SEoutput(i,ii,:) > max(max(SEoutput))*n);
%         dist1(i,ii,:) = activated;
%     end
%     
% end
% 
% for i = 1:size(SEoutput,1) % For every electrode
%     
%     for ii = 1:size(SEoutput,3) % For every neuron
%         % Find the first activated unit step
%         activated = dist1(i,:,ii);
%         if any(activated == 1)
%             dist2(i,ii) = units(find(activated == 1,1));
%         else
%             dist2(i,ii) = nan;
%         end
%     end
%     
% end
% 
% figure;
% hist(dist2(45,:),50); title('Distribution of 50% Current Intensity Threshold Values for center electrode'); xlabel('Current Intensity mA/mm^2'); ylabel('Number of Neurons Activated');

end