function Neuro = RunBaseline(Params,Neuro)
% Neuro = RunBaseline(Params,Neuro)
% Baseline period - grab neural features to get baseline for z-scoring
% 
% Neuro - contains NeuralStats structure w/ mean and var of each chan

fprintf('Collecting Baseline')

tstart  = GetSecs;
tlast = tstart;
done = 0;
while ~done,
    % Update Time & Position
    tim = GetSecs;
    
    % for pausing and quitting expt
    if CheckPause, ExperimentPause(Params); end
    
    % Grab data every Xsecs
    if (tim-tlast) > 1/Params.ScreenRefreshRate,
        % time
        tlast = tim;
        
        % grab and process neural data
        if Params.BLACKROCK && ((tim-Neuro.LastUpdateTime)>1/Params.NeuralRefreshRate),
            Neuro.LastUpdateTime = tim;
            Neuro = NeuroPipeline(Neuro);
            % update command line with progress
            fprintf('.')
        end
        
        % update screen with progress
        tex = sprintf('Computing Baseline: %.1f%% ', 100*(tim-tstart)/Params.BaselineTime);
        DrawFormattedText(Params.WPTR, tex,'center','center',255);
        Screen('Flip', Params.WPTR);
        
    end
    
    % end if takes too long
    if (tim - tstart) > Params.BaselineTime,
        done = 1;
    end
end

fprintf('Done\n\n')

end % RunBaseline