% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      222222   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          22 22    %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE         22     %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE           22      %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      222222   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;               % clear workspace and figures

d = 5;                                      % # parameters?
N = 5000;                                   % # samples to be used
mu = 0.5*ones(1,d);                         % Sample mean, µ
C = eye(d) + diag([0.6 0.2 0 0.2],-1) ...   % Covariance matrix,Σ: Eq. 44
        + diag([0.6 0.2 0 0.2],1) ... 
        + diag([0.2 0 0],2) ...
        + diag([0.2 0 0],-2);
% C = [ 1.0 0.6 0.2 0.0 0.0 ;
%       0.6 1.0 0.2 0.0 0.0 ;
%       0.2 0.2 1.0 0.0 0.0 ;
%       0.0 0.0 0.0 1.0 0.2 ;
%       0.0 0.0 0.0 0.2 1.0 ]'
switch isfile('Xy.mat')
    case 0 % Create data
        X = mvnrnd(mu,C,N);                 % Draw N samples from N(µ,Σ) 
        X = (X - repmat(min(X),N,1))./ ...  % Normalize X - values
            repmat(max(X)-min(X),N,1);              
        y = sum(X,2);                       % y = f(x) function
        save Xy.mat X y                     % Save to file
    otherwise
        load Xy.mat
end
sigma2 = var(y)/100;                        % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR_EXT options
    'basis',1,'maxorder',3,'m',3, ...
    'K',1,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Now run the HDMR toolbox
