function ExperimentStart(Subject,ControlMode,BLACKROCK,DEBUG)
% function ExperimentStart(Subject,ControlMode)
% Subject - string for the subject id
% ControlMode - [1,2,3] for mouse pos, mouse vel, & kalman
% BLACKROCK - [0,1] if 1, collects, processes, and saves neural data
% DEBUG - [0,1] if 1, enters DEBUG mode in which screen is small and cursor
%   remains unhidden

%% Clear All and Close All
clearvars -global -except Subject ControlMode BLACKROCK DEBUG
clc
warning off

if ~exist('Subject','var'), Subject = 'Test'; DEBUG = 1; end
if ~exist('ControlMode','var'), ControlMode = 2; end
if ~exist('BLACKROCK','var'), BLACKROCK = 0; end
if ~exist('DEBUG','var'), DEBUG = 0; end

AssertOpenGL;
KbName('UnifyKeyNames');

if strcmpi(Subject,'Test'), Subject = 'Test'; end

%% Retrieve Parameters from Params File
Params.Subject = Subject;
Params.ControlMode = ControlMode;
Params.BLACKROCK = BLACKROCK;
Params.DEBUG = DEBUG;
Params = GetParams(Params);

%% Initialize Blackrock System
if BLACKROCK,
    addpath('C:\Program Files (x86)\Blackrock Microsystems\NeuroPort Windows Suite')
    cbmex('close'); % always close
    cbmex('open'); % open library
    cbmex('trialconfig', 1); % empty the buffer
end

%% Neural Signal Processing
% create neuro structure for keeping track of all neuro updates/state
% changes
Neuro.ZscoreRawFlag     = Params.ZscoreRawFlag;
Neuro.ZscoreFeaturesFlag= Params.ZscoreFeaturesFlag;
Neuro.DimRed            = Params.DimRed;
Neuro.CLDA              = Params.CLDA;
Neuro.SaveProcessed     = Params.SaveProcessed;
Neuro.FilterBank        = Params.FilterBank;
Neuro.NumChannels       = Params.NumChannels;
Neuro.BufferSamps       = Params.BufferSamps;
Neuro.BadChannels       = Params.BadChannels;
Neuro.ReferenceMode     = Params.ReferenceMode;
Neuro.NumFeatures       = Params.NumFeatures;
Neuro.LastUpdateTime    = GetSecs;

% initialize filter bank state
for i=1:length(Params.FilterBank),
    Neuro.FilterBank(i).state = [];
end

% initialize stats for each channel for z-scoring
Neuro.ChStats.wSum1  = 0; % count
Neuro.ChStats.wSum2  = 0; % squared count
Neuro.ChStats.mean   = zeros(1,Params.NumChannels); % estimate of mean for each channel
Neuro.ChStats.S      = zeros(1,Params.NumChannels); % aggregate deviation from estimated mean for each channel
Neuro.ChStats.var    = zeros(1,Params.NumChannels); % estimate of variance for each channel

% initialize stats for each feature for z-scoring
Neuro.FeatureStats.wSum1  = 0; % count
Neuro.FeatureStats.wSum2  = 0; % squared count
Neuro.FeatureStats.mean   = zeros(1,Params.NumChannels); % estimate of mean for each channel
Neuro.FeatureStats.S      = zeros(1,Params.NumChannels); % aggregate deviation from estimated mean for each channel
Neuro.FeatureStats.var    = zeros(1,Params.NumChannels); % estimate of variance for each channel

% create low freq buffers
Neuro.FilterDataBuf = zeros(Neuro.BufferSamps,Neuro.NumChannels,3);

%% Kalman Filter
if Params.ControlMode==3,
    KF = Params.KF;
    KF.CLDA = Params.CLDA;
else,
    KF = [];
end

%% Check Important Params with User
LogicalStr = {'false', 'true'};
IMStr = {'imagined mvmts', 'shuffled imagined mvmts'};
DimRedStr = {'PCA', 'FA'};
Params.Subject = Subject;
Params.ControlMode = ControlMode;
Params.BLACKROCK = BLACKROCK;
Params.DEBUG = DEBUG;

fprintf('\n\nImportant Experimental Parameters:')
fprintf('\n\n  Task Parameters:')
fprintf('\n    - task: %s', Params.Task)
fprintf('\n    - subject: %s', Params.Subject)
fprintf('\n    - control mode: %s', Params.ControlModeStr)
fprintf('\n    - blackrock mode: %s', LogicalStr{Params.BLACKROCK+1})
fprintf('\n    - debug mode: %s', LogicalStr{Params.DEBUG+1})

fprintf('\n\n  Neuro Processing Pipeline:')
if Params.GenNeuralFeaturesFlag,
    fprintf('\n    - generating neural features!')
else,
    fprintf('\n    - reference mode: %s', Params.ReferenceModeStr)
    fprintf('\n    - zscore raw: %s', LogicalStr{Params.ZscoreRawFlag+1})
    fprintf('\n    - zscore features: %s', LogicalStr{Params.ZscoreFeaturesFlag+1})
    fprintf('\n    - save filtered data: %s', LogicalStr{Params.ZscoreRawFlag+1})
end
fprintf('\n    - dimensionality reduction: %s', LogicalStr{Params.DimRed.Flag+1})
if Params.DimRed.Flag,
    fprintf('\n      - method: %s', DimRedStr{Params.DimRed.Method})
end

fprintf('\n\n  BCI Parameters:')
fprintf('\n    - Imagined Movements: %s', LogicalStr{double(Params.NumImaginedBlocks>0) +1})
fprintf('\n      - initialization mode: %s', IMStr{Params.InitializationMode})
fprintf('\n    - Adaptation Decoding: %s', LogicalStr{double(Params.NumAdaptBlocks>0) +1})
fprintf('\n      - adapt type: %s', Params.CLDA.TypeStr)
fprintf('\n      - adapt change type: %s', Params.CLDA.AdaptType)
fprintf('\n    - Fixed Decoding: %s', LogicalStr{double(Params.NumFixedBlocks>0) +1})


str = input('\n\nContinue? (''n'' to quit, otherwise continue)\n' ,'s');
if strcmpi(str,'n'),
    fprintf('\n\nExperiment Ended\n\n')
    return
end

%% Initialize Window
% Screen('Preference', 'SkipSyncTests', 0);
if DEBUG
    [Params.WPTR, Params.ScreenRectangle] = Screen('OpenWindow', 0, 0, [50 50 1000 1000]);
else
    [Params.WPTR, Params.ScreenRectangle] = Screen('OpenWindow', 0, 0, [10 30 1900 1000]);
end
Params.Center = [mean(Params.ScreenRectangle([1,3])),mean(Params.ScreenRectangle([2,4]))];
if ~DEBUG, HideCursor; end

% Font
Screen('TextFont',Params.WPTR, 'Arial');
Screen('TextSize',Params.WPTR, 28);

%% Start
try
    % Baseline 
    if Params.BaselineTime>0,
        Neuro = RunBaseline(Params,Neuro);
    end
    
    % Imagined Cursor Movements Loop
    if Params.NumImaginedBlocks>0,
        [Neuro,KF] = RunTask(Params,Neuro,1,KF);
    end
    
    % Adaptation Loop
    if Params.NumAdaptBlocks>0,
        [Neuro,KF] = RunTask(Params,Neuro,2,KF);
    end
    
    % Fixed Decoder Loop
    if Params.NumFixedBlocks>0,
        [Neuro,KF] = RunTask(Params,Neuro,3,KF);
    end
    
    % Pause and Finish!
    ExperimentStop();
    
catch ME, % handle errors gracefully
    Screen('CloseAll')
    for i=length(ME.stack):-1:1,
        if i==1,
            errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s\n\n', ...
                ME.stack(1).name, ME.stack(1).line, ME.message);
        else,
            errorMessage = sprintf('Error in function %s() at line %d.\n\n', ...
                ME.stack(i).name, ME.stack(i).line);
        end
        fprintf(1,'\n%s\n', errorMessage);
    end
    keyboard;
end

end % ExperimentStart
