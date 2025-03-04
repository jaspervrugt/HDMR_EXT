function [Xy,Em,S,O,n_ns,ct] = HDMR_EXT_initialize(X,y,N,d,K,R, ...
    m,maxorder,basis)
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
% Initialize main variables/structures of HDMR_EXT                        %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Xy,Em,S,O,n_ns,ct] = HDMR_EXT_initialize(X,y,N,d,K,R, ...            %
%       m,maxorder,basis)                                                 %
%  where                                                                  %
%    X           Nxd matrix: N vectors of d parameters                    %
%    y           Nx1 vector: single model output for each row of X        %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Random seed (legacy: randn('state',sum(100*clock)); )
rng(1+round(100*rand),'twister');

% STRUCTURE XY: Define content
if K == 1
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1), ...
        'maxX',max(X,[],1),'y',y,'R',R,'id',(1:R)');
else
    % Now setup the boostrap matrix with selection matrix, id, for samples
    [~,id] = sort(rand(N,K));
    % Store in structure XY
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1), ...
        'maxX',max(X,[],1),'y',y,'R',R,'id',id(1:R,1:K));
end

% Define a and b for different basis functions
switch basis
    case 1  % orthonormal polynomial
        a = 0; b = 1;
    case 2  % legendre polynomial
        a = -1; b = 1;
    case 3  % hermite polynomial
        a = -2; b = 2;
end
% Compute normalized X-values
Xy.X_n = a + (b - a) * (X(1:N,1:d) - repmat(Xy.minX(1:d),N,1)) ./ ...
    repmat(Xy.maxX(1:d)-Xy.minX(1:d),N,1);

% Define matrix C
switch basis
    case 1  % Classical Orthonormal Polynomials
        C = orthopoly(Xy.X_n,N,d,m);
    case 2  % Legendre Polynomials
        C = legendre(m,d);
    case 3  % Hermite Polynomials
        C = hermite(m,d);
end
% Determine all combinations of parameters (coefficients) for each order
n1 = d;                             % # component functions in 1st order
[n2,n3,beta,gamma,ct] = deal(0);    % Deal zeros to n2, n3, beta, gamma and ct
n_ns = nan(1,K);                    % Pre-allocate # of nonsingular values
if maxorder > 1
    beta = nchoosek(1:d,2); n2 = size(beta,1);
end
if maxorder == 3
    gamma = nchoosek(1:d,3); n3 = size(gamma,1);
end
n = n1 + n2 + n3;                   % Total number of component functions
m1 = m;                             % # terms in 1st order for one term
m2 = 2*m + m^2;                     % # terms in 2nd order for one term
m3 = 3*m + 3*m^2 + m^3;             % # terms in 3rd order for one term
k1 = d * m;                         % Total number of terms in 1st order
k2 = m2 * (d*(d-1)/2);              % Total number of terms in 2nd order
k3 = m3 * (d*(d-1)*(d-2)/6);        % Total number of terms in 3rd order
k = k1 + k2 + k3;                   % Total number of terms

% STRUCTURE Em: Initialization
switch maxorder
    case 1
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'm1',m1,'n1',d,'m2',m2,'n2',n2,'m3',m3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'k1',k1,'k2',k2,'k3',k3,'k',k1,'beta',beta,'gamma',gamma, ...
            'C',nan(k1,K),'EB',zeros(N,k1),...
            'RT',zeros(K,1));
    case 2
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'm1',m1,'n1',n1,'m2',m2,'n2',n2,'m3',m3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'k1',k1,'k2',k2,'k3',k3,'k',k1+k2,'beta',beta,'gamma',gamma, ...
            'C',nan(k1+k2,K),'EB',zeros(N,k1+k2),...
            'RT',zeros(K,1));
    case 3
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'm1',m1,'n1',n1,'m2',m2,'n2',n2,'m3',m3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'k1',k1,'k2',k2,'k3',k3,'k',k1+k2+k3,'beta',beta, ...
            'gamma',gamma,'C',nan(k,K),'EB',zeros(N,k),...
            'RT',zeros(K,1));
end

% Now compute B-spline values for all N samples of X_n
[Em.EB,O] = HDMR_EXT_construct_EB(Em.EB,Xy.X_n(1:N,1:d),N,d,m,C,maxorder);
%Em.EB(4991:5000,end-9:end)
%pause
% STRUCTURE SA: Sensitivity analysis and analysis of variance decomposition
S = struct('S',nan(Em.n,K),'Sa',nan(Em.n,K),'Sb',nan(Em.n,K), ...
    'ST',nan(d,1),'V_em',nan(Em.n,K),'V_y',nan(1,K));

% Print to screen
fprintf('  -----------------------------------------------------------------------------  \n');
fprintf('     HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR     EEEEEEEE XXX   XXX TTTTTTTTT     \n');
fprintf('     HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT     \n');
fprintf('     HHH   HHH DDD   DDD MMM    MMM RRR   RRR   EEE       XXX XXX  TT TTT TT     \n');
fprintf('     HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR   EEE       XXX XXX  T  TTT  T     \n');
fprintf('     HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR    EEEEEE     XXXXX      TTT        \n');
fprintf('     HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR     EEEEEE     XXXXX      TTT        \n');
fprintf('     HHH   HHH DDD   DDD MMM    MMM RRRRRRR     EEE       XXX XXX     TTT        \n');
fprintf('     HHH   HHH DDD   DDD MMM    MMM RRR  RRR    EEE       XXX XXX     TTT        \n');
fprintf('     HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT        \n');
fprintf('     HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR   EEEEEEEE XXX   XXX    TTT        \n');
fprintf('  -----------------------------------------------------------------------------  \n');
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n'); 

end
