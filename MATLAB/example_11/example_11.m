% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE    11  11     %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE        11  11     %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE     11  11     %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        11  11     %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE    11  11     %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;               % clear workspace and figures

% Modified Sobol test function
d = 10;                                     % # parameters? 
N = 5000;                                   % # samples to be used
X = rand(N,d);                              % Sample X-values from U(0,1)
y = nan(N,1);                               % Initialize y values
for i = 1:N
    y(i,1) = test_function(X(i,1:d),d);     % Evaluate modified Sobol func
end

options = struct('graphics',1, ...          % Specify HDMR_EXT options
    'basis',3,'maxorder',2,'m',3, ...
    'K',10,'R',N/5,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Now run the HDMR_EXT toolbox