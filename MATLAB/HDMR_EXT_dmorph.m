function [Cf,n_ns] = HDMR_EXT_dmorph(EB,B,C,Y_res,R,k1,k)
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
% Function executes D-MORPH regression analysis                           %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Cf,n_ns] = HDMR_EXT_dmorph(EB,B,C,Y_res,R,k1,k)                      %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

C(1:k1,:) = [];                 % Remove all 1st order component functions
D = EB' * Y_res / R;
D(1:k1,:) = [];                 % Right hand size of the algebraic equation

C_inv = pinv(C);                % Generalized inverse of coefficient matrix
Cf_0 = C_inv * D;               % Least squares solution

Ia = eye(k);
Pr = Ia - C_inv * C;
PB = Pr * B;

[U,S,V] = svd(PB);              % Singular value decomposition
s = diag(S);                    % Extract vector with singular values
diff_s = -diff(s) ./ s(2:k);    % Normlzd diff adjacent singular values
tol = 0;                        % Initialize tolerance
n_ns = 1;                       % Index non-singular values not satify tol 
for i = 1 : k-1
    if diff_s(i) > tol
        tol = diff_s(i);
        n_ns = i;
    end
end

U(:,1:n_ns) = [];               % Remove columns singular values from U
V(:,1:n_ns) = [];               % Remove columns singular values from V
UV = U' * V;
Q = V * pinv(UV) * U';
Cf = Q * Cf_0;                  % D-Morph Regression

end
