
function [] = drawNoiseCorrChangeByTime(cross_noise_correlation,c,frames,color)
categorized_data  = squeeze(cross_noise_correlation(:,:,c,:,:));
reshaped_categorized_data = reshape(categorized_data,[],frames);
datamean = mean(reshaped_categorized_data,1,'omitnan');
datastd = std(reshaped_categorized_data,1,'omitnan');
datasem = datastd./sqrt(size(reshaped_categorized_data,1));

errorshade([1:frames],datamean,datasem,'lineprops',{'Color',color});
end