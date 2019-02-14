function Neuro = CompNeuralFeatures(Neuro),
% Neuro = CompNeuralFeatures(Neuro)
% computes neural features
% phase in delta band + pwr in all available bands
% sets features on bad channels to 0
%
% Neuro
%   .DeltaBuf - buffer of delta band filtered neural data [ samps x chans ]
%   .FilteredData - filtered data from last bin [ samps x chans x frqs ]
%   .NeuralFeatures - vector of features for decoding [ features*chans x 1 ]

% allocate memory
samps = Neuro.NumSamps;
neural_features = zeros(Neuro.NumFeatures,Neuro.NumChannels);

% first compute hilbert for first 3 freq bands
H = hilbert(Neuro.FilterDataBuf);

% compute phase for delta band using hilbert (only keep last bin)
ang = angle(H(:,:,1));
neural_features(1,:) = angle(sum(exp(1i*ang(end-samps+1:end,:))));

% compute pwr in first 3 frq bands based on hilbert (only keep last bin)
pwr = abs(H);
neural_features(2:4,:) = log10(squeeze(mean(pwr(end-samps+1:end,:,:))))';

% compute average pwr for all remaining frq bands in last bin
pwr = log10(mean(Neuro.FilteredData(:,:,4:end).^2, 1));

% combine feature vectors and remove singleton dimension
feature_idx = [Neuro.FilterBank.feature];
feature_idx = feature_idx(4:end); % already computed for first 3 freq bands
for i=4:Neuro.NumFeatures-1,
    idx = feature_idx == i;
    neural_features(i+1,:) = mean(pwr(:,:,idx),3);
end

% set bad channels to 0
neural_features(:,Neuro.BadChannels) = 0;

% put features in Neuro
Neuro.NeuralFeatures = neural_features(:);

end % CompNeuralFeatures

