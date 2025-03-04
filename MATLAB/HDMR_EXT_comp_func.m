function Y_em = HDMR_EXT_comp_func(EB,w,R,d,m,k1,k2,k3,k,n1,n2,n3,maxorder)
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
% Function computes the component functions based on weight vector, w     %
%                                                                         %
%  SYNOPSIS                                                               %
%   Y_em = HDMR_EXT_comp_func(EB,w,R,d,m,k1,k2,k3,k,n1,n2,n3,maxorder)    %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Make sure that w is a row vector
w = w(:); w = w';
% Pre-allocate memory 
Y_t = zeros(R,k); 
if maxorder == 1
    Y_em = zeros(R,n1);
elseif maxorder == 2
    Y_em = zeros(R,n1+n2);
else
    Y_em = zeros(R,n1+n2+n3);
end
% Compute temporary matrix
for i = 1:R
    Y_t(i,:) = EB(i,:) .* w;
end

% First order component functions
for i = 1:n1
    Y_em(1:R,i) = sum(Y_t(1:R,i*m-m+1:i*m),2);
end
% Second order component functions
if ( maxorder > 1 )
    m = k2/n2; % JAV: WHY RECOMPUTE M? 
    for i = 1:n2
        Y_em(1:R,i+n1) = sum(Y_t(:,i*m-m+k1+1:i*m+k1),2);
    end
end
% Third order component functions
if ( maxorder == 3 )
    m = k3/n3;
    for i = 1:n3
        Y_em(1:R,i+n1+n2) = sum(Y_t(:,i*m-m+k1+k2+1:i*m+k1+k2),2);
    end
end
% disp(Y_em)
% size(Y_em)

end
