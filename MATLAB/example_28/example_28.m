% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;           % clear workspace and figures

d = 2;                                  % # parameters?
N = 5000;                               % # samples to be used
X = normrnd(0,1,N,1);
X(:,2) = X(:,1).^2; 
y = X(:,1) + X(:,2); 
sigma2 = var(y)/100;                    % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);    % Add random error, y = Nx1 vector

options = struct('graphics',1, ...      % Specify HDMR_EXT options
    'basis',1,'maxorder',2,'m',3, ...
    'K',1,'R',1000,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
% Now run the HDMR_EXT toolbox
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    
