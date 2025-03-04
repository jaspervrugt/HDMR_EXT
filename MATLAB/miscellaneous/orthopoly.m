function C = orthopoly(X,N,d,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the coefficients of an orthonormal polynomial    %
% for N x d matrix X                                                      %
%                                                                         %
%  SYNOPSIS                                                               %
%   C = orthopoly(X,N,d,m)                                                %
%                                                                         %
%  REFERENCE                                                              %
%   Szegö, G (1939), Orthogonal Polynomials, American Mathematical        %
%       Society, https://books.google.com/books?id=ZOhmnsXlcY0C           %
%                                                                         %
%  © Written by Jasper A. Vrugt, Sept. 9, 2017                            %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Now compute moment matrix, m, for each dimension of X
k = 0; M = zeros(m+1,m+1);
for i = 1:d
    for j = 1:m+1
        for z = 1:m+1
            M(j,z,i) = sum(X(1:N,i).^k)/N;
            k = k + 1;
        end
        k = j;
    end
    k = 0;      

end
% Pre-allocate coefficient matrix
C = zeros(m,m+1,d);
% Calculate coefficients of orthonormal polynomial for each dimension of X
for i = 1:d
    for j = 1:m
        for k = 1:j+1
            z = 1:j+1; z(k) = [];
            det_ij = det(M(1:j,1:j,i)) * det(M(1:j+1,1:j+1,i));
%            disp(det(M(1:j, z, i))/sqrt(det_ij)),disp(i),disp(j),disp(k),disp(z)   
            C(j,j+2-k,i) = (-1)^(j+k+1) * det(M(1:j,z,i)) / sqrt(det_ij);
        end
    end
end

end

% % % Gram-Schmidt orthogonalization
% % X = X';
% % [m,n] = size(X); 
% % Q = zeros(m,n); R = zeros(n,n);
% % for j = 1:n
% %     v = X(:,j);
% %     for i = 1:j-1
% %         R(i,j) = Q(:,i)'*X(:,j);
% %         v = v - R(i,j)*Q(:,i);
% %     end
% %     R(j,j) = norm(v);
% %     Q(:,j) = v/R(j,j);
% % end
