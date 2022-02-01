function [noise_corr_matrix] = getNoiseCorrelation_singleMovie(spikeCountMatrix, sampleSizeThreshold)
% residual noise correlation (z-scored), assuming missing data.
% Ignore missing data--calculate covariance only if there is a sufficient 
% number of samples (N > threshold amount).
%
% IMPORTANT:expected input (spike count matrix): 
%
% spikeCountMatrix : [frames x units x repeats]
%

if(~exist('sampleSizeThreshold', 'var'))
    sampleSizeThreshold         = 10;
end

[num_frames, num_units, num_repeats] = size(spikeCountMatrix);

% z-score each frame's response ([unique condition x units])
z_scored_residuals  = nan(num_frames * num_repeats, num_units);
for i = 1:num_units
    
    cell_resp           = squeeze(spikeCountMatrix(:,i,:));
    frame_resp          = reshape(cell_resp, [num_frames, num_repeats]);
    unit_z_score        = [];
    
    for j = 1:num_frames
        j_frame_resp    = frame_resp(j,:);
        z_scores        = nan(size(j_frame_resp));
        if(sum(~isnan(j_frame_resp)) > sampleSizeThreshold)
            cond_mean       = nanmean(j_frame_resp);
            cond_sd         = nanstd(j_frame_resp);
            residual        = j_frame_resp - cond_mean;
            z_scores        = residual./cond_sd;
        end
        
        if(~iscolumn(z_scores))
            z_scores  = z_scores';
        end
        unit_z_score   = cat(1, unit_z_score, z_scores);
    end
    z_scored_residuals(:,i) = unit_z_score;
end

[~, noise_corr_matrix]     = cov_missing_data(z_scored_residuals, sampleSizeThreshold);
% upper_tri     = triu(noise_corr_matrix, 1); nanmean(upper_tri(:))

