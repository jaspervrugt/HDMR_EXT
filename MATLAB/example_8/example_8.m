% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      888888   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          88  88   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       888888   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          88  88   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      888888   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;               % clear workspace and figures

d = 3;                                      % # parameters? 
N = 5000;                                   % # samples to be used
mu = 0.5*ones(1,d);                         % Sample means, µ
R = [1 0.5 0.4; 0.5 1 0.7; 0.4 0.7 1];      % Sample correlation matrix
C = corr2cov([1 1 1],R);                    % Sample covariance matrix
X = mvnrnd(mu,C,N);                         % Draw N samples from N(µ,Σ)
X = (X - repmat(min(X),N,1))./ ...          % Normalize X - values
    repmat(max(X)-min(X),N,1);
y = 2*X(:,1) + X(:,2) ...                   % y = f(x) function
        + 3*exp(X(1)*X(2)*X(3));
sigma2 = var(y)/100;                        % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR_EXT options
    'basis',1,'maxorder',2,'M',3, ...
    'K',1,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Now run the HDMR_EXT toolbox
