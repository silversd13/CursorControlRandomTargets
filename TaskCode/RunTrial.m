function [Data, Neuro, KF] = RunTrial(Data,Params,Neuro,TaskFlag,KF)
% Runs a trial, saves useful data along the way
% Each trial contains the following pieces
% 1) Inter-trial interval
% 2) Get the cursor to the start target (center)
% 3) Hold position during an instructed delay period
% 4) Get the cursor to the reach target (different on each trial)
% 5) Feedback

global Cursor

%% Set up trial
ReachTargetPos = Data.TargetPosition;

% Output to Command Line
fprintf('\nTrial: %i\n',Data.Trial)
fprintf('Target: %i\n',Data.TargetPosition)
if Params.Verbose,
    fprintf('  Cursor Assistance: %.2f\n',Cursor.Assistance)
    if Params.CLDA.Type==3,
        fprintf('  Lambda: %.5g\n',Neuro.CLDA.Lambda)
    end
end

% keep track of update times
dt_vec = [];
dT_vec = [];

%% Inter Trial Interval
if ~Data.ErrorID && Params.InterTrialInterval>0,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Inter Trial Interval';

    if TaskFlag==1,
        OptimalCursorTraj = ...
            GenerateCursorTraj(Cursor.State,Cursor.State,Params.InterTrialInterval,Params);
        ct = 1;
    end
    
    done = 0;
    TotalTime = 0;
    while ~done,
        % Update Time & Position
        tim = GetSecs;

        % for pausing and quitting expt
        if CheckPause, [Neuro,Data] = ExperimentPause(Params,Neuro,Data); end

        % Update Screen Every Xsec
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate,
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt; %#ok<*AGROW>
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;
            
            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate),
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                if Params.BLACKROCK,
                    [Neuro,Data] = NeuroPipeline(Neuro,Data);
                    Data.NeuralTime(1,end+1) = tim;
                elseif Params.GenNeuralFeaturesFlag,
                    Neuro.NeuralFeatures = VelToNeuralFeatures(Params);
                    if Neuro.DimRed.Flag,
                        Neuro.NeuralFactors = Neuro.DimRed.F(Neuro.NeuralFeatures);
                        Data.NeuralFactors{end+1} = Neuro.NeuralFactors;
                    end
                    Data.NeuralFeatures{end+1} = Neuro.NeuralFeatures;
                    Data.NeuralTime(1,end+1) = tim;
                end
                KF = UpdateCursor(Params,Neuro,TaskFlag,Cursor.State(1:2),KF);
            end
            
            % cursor
            if TaskFlag==1, % imagined movements
                Cursor.State(3:4) = (OptimalCursorTraj(ct,:)'-Cursor.State(1:2))/dt;
                Cursor.State(1:2) = OptimalCursorTraj(ct,:);
                ct = ct + 1;
            end
            CursorRect = Params.CursorRect;
            CursorRect([1,3]) = CursorRect([1,3]) + Cursor.State(1) + Params.Center(1); % add x-pos
            CursorRect([2,4]) = CursorRect([2,4]) + Cursor.State(2) + Params.Center(2); % add y-pos
            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;

            % draw
            Screen('FillOval', Params.WPTR, Params.CursorColor, CursorRect);
            Screen('DrawingFinished', Params.WPTR);
            Screen('Flip', Params.WPTR);
        end

        % end if takes too long
        if TotalTime > Params.InterTrialInterval,
            done = 1;
        end

    end % Inter Trial Interval
end % only complete if no errors

%% Go to reach target
if ~Data.ErrorID,
    tstart  = GetSecs;
    Data.Events(end+1).Time = tstart;
    Data.Events(end).Str  = 'Reach Target';

    if TaskFlag==1,
        OptimalCursorTraj = [...
            GenerateCursorTraj(Cursor.State(1:2),ReachTargetPos,1.5,Params);
            GenerateCursorTraj(ReachTargetPos,ReachTargetPos,Params.TargetHoldTime,Params)];
        ct = 1;
    end
    
    done = 0;
    TotalTime = 0;
    InTargetTotalTime = 0;
    while ~done,
        % Update Time & Position
        tim = GetSecs;

        % for pausing and quitting expt
        if CheckPause, [Neuro,Data] = ExperimentPause(Params,Neuro,Data); end

        % Update Screen
        if (tim-Cursor.LastPredictTime) > 1/Params.ScreenRefreshRate,
            % time
            dt = tim - Cursor.LastPredictTime;
            TotalTime = TotalTime + dt;
            dt_vec(end+1) = dt;
            Cursor.LastPredictTime = tim;
            Data.Time(1,end+1) = tim;

            % grab and process neural data
            if ((tim-Cursor.LastUpdateTime)>1/Params.UpdateRate),
                dT = tim-Cursor.LastUpdateTime;
                dT_vec(end+1) = dT;
                Cursor.LastUpdateTime = tim;
                if Params.BLACKROCK,
                    [Neuro,Data] = NeuroPipeline(Neuro,Data);
                    Data.NeuralTime(1,end+1) = tim;
                elseif Params.GenNeuralFeaturesFlag,
                    Neuro.NeuralFeatures = VelToNeuralFeatures(Params);
                    if Neuro.DimRed.Flag,
                        Neuro.NeuralFactors = Neuro.DimRed.F(Neuro.NeuralFeatures);
                        Data.NeuralFactors{end+1} = Neuro.NeuralFactors;
                    end
                    Data.NeuralFeatures{end+1} = Neuro.NeuralFeatures;
                    Data.NeuralTime(1,end+1) = tim;
                end
                KF = UpdateCursor(Params,Neuro,TaskFlag,ReachTargetPos,KF);
            end
            
            % cursor
            if TaskFlag==1, % imagined movements
                Cursor.State(3:4) = (OptimalCursorTraj(ct,:)'-Cursor.State(1:2))/dt;
                Cursor.State(1:2) = OptimalCursorTraj(ct,:);
                ct = ct + 1;
            end
            CursorRect = Params.CursorRect;
            CursorRect([1,3]) = CursorRect([1,3]) + Cursor.State(1) + Params.Center(1); % add x-pos
            CursorRect([2,4]) = CursorRect([2,4]) + Cursor.State(2) + Params.Center(2); % add y-pos
            Data.CursorState(:,end+1) = Cursor.State;
            Data.IntendedCursorState(:,end+1) = Cursor.IntendedState;
            Data.CursorAssist(1,end+1) = Cursor.Assistance;

            % reach target
            ReachRect = Params.TargetRect; % centered at (0,0)
            ReachRect([1,3]) = ReachRect([1,3]) + ReachTargetPos(1) + Params.Center(1); % add x-pos
            ReachRect([2,4]) = ReachRect([2,4]) + ReachTargetPos(2) + Params.Center(2); % add y-pos

            % draw
            inFlag = InTarget(Cursor,ReachTargetPos,Params.TargetSize);            
            if inFlag, ReachCol = Params.InTargetColor;
            else, ReachCol = Params.OutTargetColor;
            end
            Screen('FillOval', Params.WPTR, ...
                cat(1,ReachCol,Params.CursorColor)', ...
                cat(1,ReachRect,CursorRect)')
            Screen('DrawingFinished', Params.WPTR);
            Screen('Flip', Params.WPTR);
            
            % start counting time if cursor is in target
            if inFlag,
                InTargetTotalTime = InTargetTotalTime + dt;
            else
                InTargetTotalTime = 0;
            end
        end

        % end if takes too long
        if TotalTime > Params.MaxReachTime,
            done = 1;
            Data.ErrorID = 3;
            Data.ErrorStr = 'ReachTarget';
            fprintf('\nERROR: %s\n',Data.ErrorStr)
        end

        % end if in start target for hold time
        if InTargetTotalTime > Params.TargetHoldTime,
            done = 1;
        end
    end % Reach Target Loop
end % only complete if no errors


%% Completed Trial - Give Feedback
Screen('Flip', Params.WPTR);

% output update times
if Params.Verbose,
    fprintf('Screen Update Frequency: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.ScreenRefreshRate,mean(1./dt_vec),std(1./dt_vec))
    fprintf('System Update Frequency: Goal=%iHz, Actual=%.2fHz (+/-%.2fHz)\n',...
        Params.UpdateRate,mean(1./dT_vec),std(1./dT_vec))
end

% output feedback
if Data.ErrorID==0,
    fprintf('\nSUCCESS\n')
    if Params.FeedbackSound,
        sound(Params.RewardSound,Params.RewardSoundFs)
    end
else
    if Params.FeedbackSound,
        sound(Params.ErrorSound,Params.ErrorSoundFs)
    end
    WaitSecs(Params.ErrorWaitTime);
end

end % RunTrial



