% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   11   0000   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       11  00  00  %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE    11  00  00  %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE       11  00  00  %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   11   0000   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
% Example 2 of the following paper
%  Chastaing, G., F. Gamboa, and C. Prieur (2012), Generalized Hoeffding-
%      Sobol decomposition for dependent variables - application to 
%      sensitivity analysis, Electron. J. Statistics, 6, pp. 2420-2448, 
%      https://doi.org/10.1214/12-EJS749

clc; clear; close all hidden;               % clear workspace and figures

d = 4;                                      % # parameters? 
N = 1000;                                   % # samples to be used
[mu1,mu2] = deal(zeros(2,1));               % sample means µ1 & µ2 normal mixture
w1 = 0.2;                                   % weight of first normal
r12_1 = 0.4; r12_2 = 0.37;                  % Covariance X1,X2 1st&2nd normal
C1 = [ 0.5 , r12_1 ; r12_1 , 0.5 ];         % Covariance matrix first normal
C2 = [ 0.7 , r12_2 ; r12_2 , 0.3 ];         % Covariance matrix 2nd normal
X = user_draw2(mu1,mu2,C1,C2,w1,N);         % Draw N samples X normal mixture
y = 5*X(:,1) + 4*X(:,2) + 3*X(:,3) ...      % y = f(X) function
        + 2*X(:,4);     
X = (X - repmat(min(X),N,1))./ ...          % Normalize X - values
    repmat(max(X)-min(X),N,1);
sigma2 = var(y)/100;                        % Variance of random error  
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR_EXT options
    'basis',3,'maxorder',2,'m',3, ...
    'K',50,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Now run the HDMR_EXT toolbox
