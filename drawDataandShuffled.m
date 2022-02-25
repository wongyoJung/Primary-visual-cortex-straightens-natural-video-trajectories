function [] = drawDataandShuffled(angles,shuffled_angles,frames)
% natural images
natural_frames_mean_arr = [];
natural_frames_sem_arr = [];
for f=1:frames
    
    reshaped_angles = reshape(angles(:,:,1,:,f),[],1);
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    natural_frames_mean_arr = [natural_frames_mean_arr; natural_frames_mean];
    natural_frames_sem_arr = [natural_frames_sem_arr; natural_frames_sem];

end

% unnatural images
unnatural_frames_mean_arr = [];
unnatural_frames_sem_arr = [];

for f=1:frames
    
    reshaped_angles = reshape(angles(:,:,2,:,f),[],1);
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    unnatural_frames_mean_arr = [unnatural_frames_mean_arr; natural_frames_mean];
    unnatural_frames_sem_arr = [unnatural_frames_sem_arr; natural_frames_sem];
end
figure;
errorshade([1:frames],natural_frames_mean_arr,natural_frames_sem_arr,'lineprops',{'Color',[0.2 0.2 0.8]});
errorshade([1:frames],unnatural_frames_mean_arr,unnatural_frames_sem_arr,'lineprops',{'Color',[0.8 0.2 0.2]});
legend('natural','unnatural');
title('data');


xlabel('frames');
ylabel('cosΘ');


% shuffled
% natural_frames_mean_arr_shuffle = [];
% natural_frames_sem_arr_shuffle = [];
% 
% for f=1:frames
% 
%     shuffle_reshaped_angles = reshape(shuffled_angles(:,:,1,:,f),[],1);
%     natural_frames_mean = mean(shuffle_reshaped_angles,'omitnan');
%     natural_frames_std = std(shuffle_reshaped_angles,'omitnan');
%     natural_frames_sem = natural_frames_std / sqrt(length(shuffle_reshaped_angles));
%     
%     natural_frames_mean_arr_shuffle = [natural_frames_mean_arr_shuffle; natural_frames_mean];
%     natural_frames_sem_arr_shuffle = [natural_frames_sem_arr_shuffle; natural_frames_sem];
% 
% end
% 
% % unnatural images
% unnatural_frames_mean_arr_shuffle = [];
% unnatural_frames_sem_arr_shuffle = [];
% for f=1:frames
%     reshaped_angles = reshape(shuffled_angles(:,:,2,:,f),[],1);
%     natural_frames_mean = mean(reshaped_angles,'omitnan');
%     natural_frames_std = std(reshaped_angles,'omitnan');
%     natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
%     
%     unnatural_frames_mean_arr_shuffle = [unnatural_frames_mean_arr_shuffle; natural_frames_mean];
%     unnatural_frames_sem_arr_shuffle = [unnatural_frames_sem_arr_shuffle; natural_frames_sem];
% end
% figure;
% errorshade([1:frames],natural_frames_mean_arr_shuffle,natural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.2 0.2 0.8]});
% errorshade([1:frames],unnatural_frames_mean_arr_shuffle,unnatural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.8 0.2 0.2]});
% legend('natural','unnatural');
% title('shuffled');
% xlabel('frames');
% ylabel('cosΘ');
end
