function [N,d,graphics,basis,maxorder,m,K,R,method,alfa] = ...
    HDMR_EXT_setup(X,y,user_options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT %
%  HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR   EEEEEEEE XXX   XXX TTTTTTTTT %
%  HHH   HHH DDD   DDD MMM    MMM RRR   RRR  EEE       XXX XXX  TT TTT TT %
%  HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR  EEE       XXX XXX  T  TTT  T %
%  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR   EEEEEE     XXXXX      TTT    %
%  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR    EEEEEE     XXXXX      TTT    %
%  HHH   HHH DDD   DDD MMM    MMM RRRRRRR    EEE       XXX XXX     TTT    %
%  HHH   HHH DDD   DDD MMM    MMM RRR  RRR   EEE       XXX XXX     TTT    %
%  HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    %
%  HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Check/setup of the input arguments of HDMR_EXT function                 %
%                                                                         %
%  SYNOPSIS                                                               %
%   [N,d,graphics,basis,maxorder,m,K,R,method,alfa] = ...                 %
%       HDMR_EXT_setup(X,y,options)                                       %
%  where                                                                  %
%    X           Nxd matrix: N vectors of d parameters                    %
%    y           Nx1 vector: single model output for each row of X        %
%    options     (optional) structure: Fields specify HDMR variables      %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Determine the size of matrix X
[N,d] = size(X); 
if N < 300
    error(['HDMR_EXT ERROR: Number of samples, N, of N x d matrix ' ...
        'X is insufficient\n']); 
end
[a,b] = size(y);
if a ~= N
   error(['HDMR_EXT ERROR: Dimension mismatch. The number of rows, N ' ...
       'of Y should match number of rows of X\n']); 
end 
if b ~= 1
   error(['HDMR_EXT ERROR: Y should be a N x 1 vector with one ' ...
       'simulated output for each N parameter vectors\n']); 
end 

% Define default options
def_options = struct('graphics',1,'basis',1,'maxorder',2,'m',3, ...
    'K',1,'R',floor(N/2),'alfa',0.01,'method',1);

% Open file for settings.txt
fid = fopen('HDMR_EXT_settings.txt','w');
fprintf(fid,'|-------------------------------------------------------------------------|\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'| HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR     EEEEEEEE XXX   XXX TTTTTTTTT |\n');
fprintf(fid,'| HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT |\n');
fprintf(fid,'| HHH   HHH DDD   DDD MMM    MMM RRR   RRR   EEE       XXX XXX  TT TTT TT |\n');
fprintf(fid,'| HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR   EEE       XXX XXX  T  TTT  T |\n');
fprintf(fid,'| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR    EEEEEE     XXXXX      TTT    |\n');
fprintf(fid,'| HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR     EEEEEE     XXXXX      TTT    |\n');
fprintf(fid,'| HHH   HHH DDD   DDD MMM    MMM RRRRRRR     EEE       XXX XXX     TTT    |\n');
fprintf(fid,'| HHH   HHH DDD   DDD MMM    MMM RRR  RRR    EEE       XXX XXX     TTT    |\n');
fprintf(fid,'| HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n');
fprintf(fid,'| HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT    |\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|  SSSSSSS EEEEEEEE TTTTTTTT TTTTTTTT IIIIIIII NN    NN  GGGGGG   SSSSSSS |\n');
fprintf(fid,'| SSSSSSS  EEEEEEEE TTTTTTTT TTTTTTTT  IIIIII  NNN   NN GGGGGG   SSSSSSS  |\n');
fprintf(fid,'| SS       EE       TT TT TT TT TT TT    II    NNNN  NN GG       SS       |\n');
fprintf(fid,'| SSSSSS   EEEEE    T  TT  T T  TT  T    II    NN NN NN GG  GGGG SSSSSS   |\n');
fprintf(fid,'| SSSSSSS  EEEEE       TT       TT       II    NN NN NN GG   GGG SSSSSSS  |\n');
fprintf(fid,'|       SS EE          TT       TT       II    NN  NNNN GG    GG       SS |\n');
fprintf(fid,'|  SSSSSSS EEEEEEEE    TT       TT     IIIIII  NN   NNN GGGGGGGG  SSSSSSS |\n');
fprintf(fid,'| SSSSSSS  EEEEEEEE    TT       TT    IIIIIIII NN    NN  GGGGGG  SSSSSSS  |\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|-------------------------------------------------------------------------|\n');

% Now unpack structure user_options
if isfield(user_options,'graphics')
    graphics = user_options.graphics;
    if ~any(ismember(graphics,[0 1]))
        error(['HDMR_EXT ERROR: Field "graphics" of options should ' ...
            'take on the value of 0 or 1']);
    end
else
    fprintf(['HDMR_EXT DEFAULT: Field "graphics" of options not ' ...
        'specified. Default: graphics = 1 assumed\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "graphics" of options not ' ...
        'specified. Default: graphics = 1 assumed\n']);
    graphics = def_options.graphics;
end

if isfield(user_options,'maxorder')
    maxorder = user_options.maxorder;
    if ~any(ismember(maxorder,[1 2 3]))
        error(['HDMR_EXT ERROR: Field "maxorder" of options should ' ...
            'be an integer with values of 1, 2 or 3']);
    end
else
    fprintf(['HDMR_EXT DEFAULT: Field "maxorder" of options not ' ...
        'specified. Default: maxorder = 3 used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "maxorder" of options not ' ...
        'specified. Default: maxorder = 3 used\n']);
    maxorder = def_options.maxorder;
end
if isfield(user_options,'basis')
    basis = user_options.basis;
    if ~any(ismember(basis,[1 2 3]))
        error(['HDMR_EXT ERROR: Field "basis" of options should be an ' ...
            'integer with values of 1, 2 or 3']);
    end
else
    fprintf(['HDMR_EXT DEFAULT: Field "basis" of options not ' ...
        'specified. Default: basis = 1 used {Classical Orthonormal ' ...
        'Polynomials}\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "basis" of options not ' ...
        'specified. Default: basis = 1 used {Classical Orthonormal ' ...
        'Polynomials}\n']);
    basis = def_options.basis;
end
if isfield(user_options,'m')
    m = user_options.m;
    if ~any(ismember(m,1:12))
        error(['HDMR_EXT ERROR: Field "m" of options should ' ...
            'be an integer between 1 to 12']);
    end    
else
    fprintf(['HDMR_EXT DEFAULT: Field "m" of options not ' ...
        'specified. Default: m = 3 is used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "m" of options not ' ...
        'specified. Default: m = 3 is used\n']);
    m = def_options.m;
end
if isfield(user_options,'K')
    K = user_options.K;
    if ~any(ismember(K,1:500))
        error(['HDMR_EXT ERROR: Field "K" of options should be ' ...
            'an integer between 1 to 500']);
    end        
else
    fprintf(['HDMR_EXT DEFAULT: Field "K" of options not ' ...
        'specified. Default: K = 1 is used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "K" of options not ' ...
        'specified. Default: K = 1 is used\n']);
    K = def_options.K;
end
if isfield(user_options,'R')
    R = user_options.R;
    if ~any(ismember(R,100:N))
        error(['HDMR_EXT ERROR: Field "R" of options should be ' ...
            'an integer between 100 and N, the number of rows matrix X']);
    end            
else
    fprintf(['HDMR_EXT DEFAULT: Field "R" of options not ' ...
        'specified. Default: R = floor(N/2) is used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "R" of options not ' ...
        'specified. Default: R = floor(N/2) is used\n']);
    R = def_options.R;
end
if isfield(user_options,'method')
    method = user_options.method;
    if ~any(ismember(method,[1 2]))
        error(['HDMR_EXT ERROR: Field "method" of options should ' ...
            'take on the value of 1 (forward selection) ' ...
            'or 2 (backward elimination)']);
    end    
else
    fprintf(['HDMR_EXT DEFAULT: Field "method" of options not ' ...
        'specified. Default: method = 1 (forward selection) is used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "method" of options not ' ...
        'specified. Default: method = 1 (forward selection) is used\n']);
    method = def_options.method;
end
if isfield(user_options,'alfa')
    alfa = user_options.alfa;
    if ischar(alfa)
        error(['HDMR_EXT ERROR: Field "alfa" (significance level) ' ...
            'of options should not be a string but a numerical value']);
    end 
    if alfa > 0.1
        error(['HDMR_EXT ERROR: Field "alfa" (significance level) ' ...
            'of options should not exceed 0.1. Default: 0.01']);
    end    
    if alfa < eps
        error(['HDMR_EXT ERROR: Field "alfa" (significance level) ' ...
            'of options should at least be equal to eps. Default: 0.01']);
    end    
else
    fprintf(['HDMR_EXT DEFAULT: Field "alfa" of options not ' ...
        'specified. Default: alfa = 0.01 is used\n']);
    fprintf(fid,['HDMR_EXT DEFAULT: Field "alfa" of options not ' ...
        'specified. Default: alfa = 0.01 is used\n']);
    alfa = def_options.alfa;
end

fprintf(fid,'\n');
fprintf(fid,'          |===================================|\n');
fprintf(fid,'          |  field of options      value      |\n');
fprintf(fid,'          |-----------------------------------|\n');
fprintf(fid,'          |    %-12s \t %8d     |\n','graphics   ',graphics);
fprintf(fid,'          |    %-12s \t %8d     |\n','basis      ',basis);
fprintf(fid,'          |    %-12s \t %8d     |\n','maxorder   ',maxorder);
fprintf(fid,'          |    %-12s \t %8d     |\n','m          ',m);
fprintf(fid,'          |    %-12s \t %8d     |\n','K          ',K);
fprintf(fid,'          |    %-12s \t %8d     |\n','R          ',R);
fprintf(fid,'          |    %-12s \t %8f     |\n','alfa       ',alfa);
fprintf(fid,'          |    %-12s \t %8d     |\n','method     ',method);
fprintf(fid,'          |===================================|\n');
fprintf(fid,'\n');
% Now close file
fclose(fid);

end
