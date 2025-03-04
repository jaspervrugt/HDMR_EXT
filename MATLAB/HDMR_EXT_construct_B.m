function [B,C] = HDMR_EXT_construct_B(EB,R,d,m,m2,m3,k1,k12,k,maxorder)
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
% This function constructs the cost matrix, B                             %
%                                                                         %
%  SYNOPSIS                                                               %
%   [B,C] = HDMR_EXT_construct_B(EB,R,d,m,m2,m3,k1,k12,k,maxorder)        %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Preallocate memory for C2, C3, C_22 and C_33
C_2 = nan(2*m+1,m2); C_22 = nan(m2,m2,d*(d-1)/2);
C_3 = nan(3*m+3*m^2,m3); C_33 = nan(m3,m3,d*(d-1)*(d-2)/6);
% Compute C_0 and C
C_0 = sum(EB)/R;
C = EB'*EB/R;

l = 0;
for i=1:d-1
    for j=i+1:d
        l=l+1;
        C_2(1,1:m2) = C_0(1,k1+(l-1)*m2+1:k1+l*m2);
        for z = 1 : 2 * m
            C_2(z+1,1:m2) = C(k1+(l-1)*m2+z,k1+(l-1)*m2+1:k1+l*m2);
        end
        C_22(1:m2,1:m2,l) = C_2'*C_2;
    end
end

if maxorder == 3                      
    l=0;
    for i=1:d-2
        for j=i+1:d-1
            for z=j+1:d
                l=l+1;
                C_3(1,1:m3) = C_0(1,k12+(l-1)*m3+1:k12+l*m3);
                for ii = 1 : 3*m + 3*m^2
%                    disp(C(k12+(l-1)*m3+ii,k12+(l-1)*m3+1:k12+l*m3))
                    C_3(ii+1,1:m3)=C(k12+(l-1)*m3+ii, ...
                        k12+(l-1)*m3+1:k12+l*m3);
                end
                C_33(1:m3,1:m3,l)=C_3'*C_3;
            end
        end
    end
else
    C_33 = 0;
end

B = zeros(k,k);

for l = 1 : d*(d-1)/2
    B(k1+(l-1)*m2+1:k1+l*m2,k1+(l-1)*m2+1:k1+l*m2) = C_22(:,:,l);
end

if ( maxorder == 3 )
    for l = 1 : d*(d-1)*(d-2)/6
        B(k12+(l-1)*m3+1:k12+l*m3,k12+(l-1)*m3+1:k12+l*m3) = C_33(:,:,l);
    end
end

end
