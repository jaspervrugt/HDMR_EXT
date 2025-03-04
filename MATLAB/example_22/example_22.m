% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example of the following paper
%  Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity 
%      Analysis With Dependent Variables: Covariance-Based Decomposition 
%      of Hydrologic Models, Water Resources Research, 59 (4), 
%      e2022WR0328346, https://doi.org/10.1029/2022WR032834

clc; clear; close all hidden;   % clear workspace and figures

d = 7;                          % Dimensionality of model
N = 1000;                       % Number of samples used by HDMR

% Parname:      Imax Smax Qsmax  alE alF Kfast Kslow
Par_info.min = [ 0.5  10    0   1e-6 -10  0     0  ];   % minimum values
Par_info.max = [ 10  1000  100  100   10  10   150 ];   % maximum values
label_par = {'$I_{\rm max}$','$S_{\rm max}$','$Q_{\rm s,max}$', ...
    '$\alpha_{\rm E}$','$\alpha_{\rm F}$','$K_{\rm f}$','$K_{\rm s}$'};

X = LH_sampling(Par_info.min, ...       % Draw N par. smples Latin hypercbe
    Par_info.max,N);    
Y(:,1) = rainfall_runoff(X(1,1:d));     % Run model once to determine n
ny = numel(Y); Y(:,2:N) = nan(ny,N-1);  % Initialize remainder Y matrix
switch license('checkout', ...
        'Distrib_Computing_Toolbox')
    case 0 % Sequential evaluation
        for i = 2:N
            Y(:,i) = rainfall_runoff( ...
                X(i,1:d));
        end
    case 1 % Parallel evaluation
        parfor i = 2:N
            x = X(i,1:d);               %#ok single parameter vector               
            Y(:,i) = ...
                rainfall_runoff(x);
        end
end

options = struct('graphics',1, ...	    % HDMR_EXT structure options
    'basis',1,'maxorder',2,'m',3, ...
    'K',1,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
for t = 1:1 %ny
    [S,Ss,Fx,Em,Xy] = ...		        % Run HDMR_EXT toolbox
        HDMR_EXT(X,Y(t,1:N)',options);
end
