% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      777777   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE              77   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE          77    %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE            77     %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE       77      %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;               % clear workspace and figures

d = 10;                                     % # parameters? 
N = 5000;                                   % # samples to be used
X = rand(N,d);                              % Sample X-values from U(0,1)
a = 100 * rand(1,d);                        % a-values of Sobol's g-function
                                            % a = [0 1 4.5 9 99 99 99 99 99 99]; 
y = sum((abs(4*X - 2) + repmat(a,N,1)) ...  % Compute y-values: Sobol g-function
    ./ (repmat(a,N,1) + 1),2);  
D_i = 1./(3*(1+a).^2); D = sum(D_i);        % Analytic variance-based estimates 
ST_an = 1/D * D_i';

options = struct('graphics',1, ...          % Specify HDMR_EXT options
    'basis',1,'maxorder',2,'m',3, ...
    'K',1,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Now run the HDMR_EXT toolbox

% Do numerical & analytic sensitivty estimates match?
err = sum(abs(str2double(S(2:d+1,7)) - ST_an))/d