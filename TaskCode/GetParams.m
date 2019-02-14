function Params = GetParams(Params)
% Experimental Parameters
% These parameters are meant to be changed as necessary (day-to-day,
% subject-to-subject, experiment-to-experiment)
% The parameters are all saved in 'Params.mat' for each experiment

%% Verbosity
Params.Verbose = true;

%% Experiment
Params.Task = 'RandomTargets';
switch Params.ControlMode,
    case 1, Params.ControlModeStr = 'MousePosition';
    case 2, Params.ControlModeStr = 'MouseVelocity';
    case 3, Params.ControlModeStr = 'KalmanVelocity';
end

%% Control
Params.Gain             = 1;
Params.CenterReset      = false;
Params.Assistance       = .1; % value btw 0 and 1, 1 full assist
Params.CLDA.Type        = 0; % 0-none, 1-refit, 2-smooth batch, 3-RML
Params.CLDA.AdaptType   = 'linear'; % {'none','linear'}, affects assistance & lambda for rml
Params.InitializationMode = 1; % 1-imagined mvmts, 2-shuffled imagined mvmts

%% Current Date and Time
% get today's date
now = datetime;
Params.YYYYMMDD = sprintf('%i',yyyymmdd(now));
Params.HHMMSS = sprintf('%02i%02i%02i',now.Hour,now.Minute,round(now.Second));

%% Data Saving

% if Subject is 'Test' or 'test' then can write over previous test
if strcmpi(Params.Subject,'Test'),
    Params.YYYYMMDD = 'YYYYMMDD';
    Params.HHMMSS = 'HHMMSS';
end

if IsWin,
    projectdir = 'C:\Users\ganguly-lab2\Documents\MATLAB\CursorControlRandomTargets';
elseif IsOSX,
    projectdir = '/Users/daniel/Projects/CursorControlRandomTargets/';
else,
    projectdir = '/home/dsilver/Projects/CursorControlRandomTargets/';
end
addpath(genpath(fullfile(projectdir,'TaskCode')));

% create folders for saving
datadir = fullfile(projectdir,'Data',Params.Subject,Params.YYYYMMDD,Params.HHMMSS);
Params.Datadir = datadir;
% if folders already exist, warn user before continuing (unless,
% subject='Test')
if exist(Params.Datadir,'dir'),
    if ~strcmpi(Params.Subject,'Test'),
        str = input([...
            '\n\nData directory already exists.',...
            '\nAre you sure you want to continue? Y/N ',...
            '\nThis directory will be overwritten if you continue. '],'s');
    else, str = 'Y';
    end
%     if strcmpi(str,'Y'), % delete all files in directory
%         recycle('on')
%         delete(fullfile(datadir,'*'))
%         delete(fullfile(datadir,'Imagined','*'))
%         delete(fullfile(datadir,'BCI_CLDA','*'))
%         delete(fullfile(datadir,'BCI_Fixed','*'))
%     end
end
mkdir(datadir);


%% Timing
Params.ScreenRefreshRate = 10; % Hz
Params.UpdateRate = 10; % Hz
Params.BaselineTime = 0; % secs

%% Targets
Params.TargetSize = 30;
Params.TargetRect = ...
    [-Params.TargetSize -Params.TargetSize +Params.TargetSize +Params.TargetSize];
Params.OutTargetColor = [0,255,0];
Params.InTargetColor = [255,0,0];
Params.Workspace = [-400,-400;400,400]; % [x0,y0;x1,y1]
Params.NewTargetDist = 100;
Params.TargetSelectionFlag  = 1; % 1-uniform from workspace, 2-1, but must be dist away from cursor
switch Params.TargetSelectionFlag,
    case 1, Params.TargetFunc = @() URand(Params.Workspace);
    case 2, Params.TargetFunc = @() URandDist(Params.Workspace,Params.NewTargetDist);
end

%% Cursor
Params.CursorColor = [0,0,255];
Params.CursorSize = 5;
Params.CursorRect = [-Params.CursorSize -Params.CursorSize ...
    +Params.CursorSize +Params.CursorSize];

%% Kalman Filter Properties
dt = 1/Params.UpdateRate;
if Params.ControlMode==3,
    Params.KF.A = [...
        1       0       dt      0       0;
        0       1       0       dt      0;
        0       0       .8      0       0;
        0       0       0       .8      0;
        0       0       0       0       1];
    Params.KF.W = [...
        0       0       0       0       0;
        0       0       0       0       0;
        0       0       500     0       0;
        0       0       0       500     0;
        0       0       0       0       0];
    Params.KF.P = eye(5);
    Params.KF.InitializationMode = Params.InitializationMode; % 1-imagined mvmts, 2-shuffled
end

%% Trial and Block Types
Params.NumImaginedBlocks    = 1;
Params.NumAdaptBlocks       = 1;
Params.NumFixedBlocks       = 1;
Params.NumTrialsPerBlock    = 10;

%% CLDA Parameters
TypeStrs                = {'none','refit','smooth_batch','rml'};
Params.CLDA.TypeStr     = TypeStrs{Params.CLDA.Type+1};

Params.CLDA.UpdateTime = 80; % secs, for smooth batch
Params.CLDA.Alpha = exp(log(.5) / (120/Params.CLDA.UpdateTime)); % for smooth batch
Params.CLDA.Lambda = exp(log(.5) / (30*Params.UpdateRate)); % for RML

switch Params.CLDA.AdaptType,
    case 'none',
        Params.CLDA.DeltaLambda = 0;
        Params.CLDA.DeltaAssistance = 0;
    case 'linear',
        FinalLambda = exp(log(.5) / (500*Params.UpdateRate));
        DeltaLambda = (FinalLambda - Params.CLDA.Lambda) ...
            / (Params.NumAdaptBlocks...
            *Params.NumTrialsPerBlock...
            *Params.UpdateRate...
            *5); % sec/trial;
        Params.CLDA.DeltaLambda = DeltaLambda; % for RML
        switch Params.CLDA.Type,
            case 2, % smooth batch
                Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                Params.Assistance...
                /(Params.NumAdaptBlocks*Params.NumTrialsPerBlock*5/Params.CLDA.UpdateTime);
            case 3, % RML
            Params.CLDA.DeltaAssistance = ... % linearly decrease assistance
                Params.Assistance...
                /(Params.NumAdaptBlocks*Params.NumTrialsPerBlock);
            otherwise, % none or refit
            Params.CLDA.DeltaAssistance = 0;
        end
end

%% Hold Times
Params.TargetHoldTime = .3;
Params.InterTrialInterval = 0;
Params.MaxStartTime = 10;
Params.MaxReachTime = 10;
Params.InterBlockInterval = 0;

%% Feedback
Params.FeedbackSound = false;
Params.ErrorWaitTime = 2;
Params.ErrorSound = 1000*audioread('buzz.wav');
Params.ErrorSoundFs = 8192;
[Params.RewardSound,Params.RewardSoundFs] = audioread('reward1.wav');
% play sounds silently once so Matlab gets used to it
sound(0*Params.ErrorSound,Params.ErrorSoundFs)

%% BlackRock Params
Params.GenNeuralFeaturesFlag = true;
Params.ZscoreRawFlag = true;
Params.ZscoreFeaturesFlag = false;
Params.SaveProcessed = false;

Params.DimRed.Flag = false;
Params.DimRed.Method = 1; % 1-pca, 2-fa
Params.DimRed.AvgTrialsFlag = false; % 0-cat imagined mvmts, 1-avg imagined mvmts
Params.DimRed.NumDims = [];

Params.Fs = 1000;
Params.NumChannels = 128;
Params.BufferTime = 2; % secs longer for better phase estimation of low frqs
Params.BufferSamps = Params.BufferTime * Params.Fs;
Params.BadChannels = [];
RefModeStr = {'none','common_mean','common_median'};
Params.ReferenceMode = 0; % 0-no ref, 1-common mean, 2-common median
Params.ReferenceModeStr = RefModeStr{Params.ReferenceMode+1};

% filter bank - each element is a filter bank
% fpass - bandpass cutoff freqs
% feature - # of feature (can have multiple filters for a single feature
% eg., high gamma is composed of multiple freqs)
Params.FilterBank = [];
Params.FilterBank(end+1).fpass = [.5,4];    % delta
Params.FilterBank(end).feature = 1;
Params.FilterBank(end+1).fpass = [4,8];     % theta
Params.FilterBank(end).feature = 2;
Params.FilterBank(end+1).fpass = [8,13];    % alpha
Params.FilterBank(end).feature = 3;
Params.FilterBank(end+1).fpass = [13,19];   % beta1
Params.FilterBank(end).feature = 4;
Params.FilterBank(end+1).fpass = [19,30];   % beta2
Params.FilterBank(end).feature = 4;
Params.FilterBank(end+1).fpass = [30,36];   % low gamma1 
Params.FilterBank(end).feature = 5;
Params.FilterBank(end+1).fpass = [36,42];   % low gamma2 
Params.FilterBank(end).feature = 5;
Params.FilterBank(end+1).fpass = [42,50];   % low gamma3
Params.FilterBank(end).feature = 5;
Params.FilterBank(end+1).fpass = [70,77];   % high gamma1
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
Params.FilterBank(end).feature = 6;
Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
Params.FilterBank(end).feature = 6;
% compute filter coefficients
for i=1:length(Params.FilterBank),
    [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
    Params.FilterBank(i).b = b;
    Params.FilterBank(i).a = a;
end

Params.NumFeatures = length(unique([Params.FilterBank.feature])) + 1;

%% Save Parameters
save(fullfile(Params.Datadir,'Params.mat'),'Params');

end % GetParams

