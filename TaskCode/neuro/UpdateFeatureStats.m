function Neuro = UpdateFeatureStats(Neuro)
% function Neuro = UpdateFeatureStats(Neuro)
% update estimate of mean and variance for each channel using Welford Alg
% Neuro 
% 	.NeuralFeatures - [ features x 1 ]
%   .FeatureStats - structure, which is updated

X = Neuro.NeuralFeatures';

% updates
w                           = 1;
Neuro.FeatureStats.wSum1    = Neuro.FeatureStats.wSum1 + w;
Neuro.FeatureStats.wSum2    = Neuro.FeatureStats.wSum2 + w*w;
meanOld                     = Neuro.FeatureStats.mean;
Neuro.FeatureStats.mean     = meanOld + (w / Neuro.FeatureStats.wSum1) * mean(X - repmat(meanOld,w,1));
Neuro.FeatureStats.S        = Neuro.FeatureStats.S + w*mean( (X - repmat(meanOld,w,1)).*(X - repmat(Neuro.FeatureStats.mean,w,1)) );
Neuro.FeatureStats.var      = Neuro.FeatureStats.S / (Neuro.FeatureStats.wSum1 - 1);

end % UpdateNeuralStats