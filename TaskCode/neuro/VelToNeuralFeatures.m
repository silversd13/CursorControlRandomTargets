function Z = VelToNeuralFeatures(Params,noise,PLOT)
% function Z = VelToNeuralFeatures(Params,noise,PLOT)
% Use a 2D gaussian function to generate neural features vector
% neural features change depending input velocity
%
% INPUT:
% noise - % noise in % of centroid peak value (default=100)
% PLOT - 0-no plot, 1-plot (default=1)
%
% OUTPUT: 
% Z - neural features vector with size of 128*7 
%
% CREATED: G. Nootz  May 2012
% 
%  Modifications:
%  Reza Abiri Feb 2019
%  Daniel Silversmith Feb 2019
% 
% ---------User Input---------------------

% inputs
if ~exist('noise','var'), noise=100; end
if ~exist('PLOT','var'), PLOT = 0; end

% compute velocities
[x,y] = GetMouse();
Vx = Params.Gain * (x - Params.Center(1));
Vy = Params.Gain * (y - Params.Center(2));

% rescaling to matrix map
MdataSizeY=32;
MdataSizeX=28;

Vy=MdataSizeY/2 + Vy*(MdataSizeY/(2*600));
Vx=MdataSizeX/2 + Vx*(MdataSizeX/(2*600));

% generate centroid
[X,Y] = meshgrid(1:MdataSizeX,1:MdataSizeY);

xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;

x = [2,Vx,7,Vy,4.5,+0.02*2*pi]; % centroid parameters
Z = D2GaussFunction(x,xdata);

% add noise
noise = noise/100 * x(1);
Z = Z + noise*(rand(size(X,1),size(Y,2))-0.5);

% feature plot
if PLOT,
	imagesc(X(1,:),Y(:,1)',Z)
	set(gca,'YDir','reverse')
	colormap('jet')
end

% vectorize output
Z = Z(:);

end % VelToNeuralFeatures

function F = D2GaussFunction(x,xdata)
F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) ...
    + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2)));
end % D2GaussFunction
 