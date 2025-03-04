function [EB,O] = HDMR_EXT_construct_EB(EB,X,N,d,m,C,maxorder)
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
% This function constructs the extended bases matrix, EB                  %
%                                                                         %
%  SYNOPSIS                                                               %
%   [EB,O] = HDMR_EXT_construct_EB(EB,X,N,d,m,C,maxorder)                 %
%                                                                         %
%  © Written by Jasper A. Vrugt & A. Sahin, 2017                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

k = 0;
for i = 1:d
    for j = 1:m
        k = k + 1;
        EB(:,k) = orthopolyeval(C(j,:,i),X(1:N,i),j);
        O(k,1) = j;
    end
end

if ( maxorder > 1 )
    for i = 1:d-1
        for j = i+1:d
            for z = 1:m
                k = k+1;
                EB(:,k) = orthopolyeval(C(z,:,i),X(1:N,i),z);
                O(k,1) = z;
            end
            for z = 1:m
                k = k+1;
                EB(:,k) = orthopolyeval(C(z,:,j),X(1:N,j),z);
                O(k,1) = z;
            end
            for z = 1:m
                for r = 1:m
                    k = k+1;
                    EB(:,k) = orthopolyeval(C(z,:,i),X(1:N,i),z) .* orthopolyeval(C(r,:,j),X(1:N,j),r);
                    O(k,1) = z; O(k,2) = r;
                end
            end
        end
    end
end

if ( maxorder == 3 )
    for i = 1:d-2
        for j = i+1:d-1
            for z = j+1:d
                for r = 1:m
                    k = k+1;
                    EB(:,k) = orthopolyeval(C(r,:,i),X(1:N,i),r);
%                    disp(EB(1:10,k)); disp(C(r,:,i)); disp(X(1:10,i)); disp(r);
                    O(k,1) = r;
                end
                for r = 1:m
                    k = k+1;
                    EB(:,k) = orthopolyeval(C(r,:,j),X(1:N,j),r);
                    O(k,1) = r;
                end
                for r = 1:m
                    k = k+1;
                    EB(:,k) = orthopolyeval(C(r,:,z),X(1:N,z),r);
                    O(k,1) = r;
                end
                for t = 1:m
                    for r = 1:m
                        k = k+1;
                        EB(:,k) = orthopolyeval(C(t,:,i),X(1:N,i),t) .* ...
                            orthopolyeval(C(r,:,j),X(1:N,j),r);
                        O(k,1) = t; O(k,2) = r;
                    end
                end
                for t = 1:m
                    for r = 1:m
                        k = k+1;
                        EB(:,k) = orthopolyeval(C(t,:,i),X(1:N,i),t) .* ...
                            orthopolyeval(C(r,:,z),X(1:N,z),r);
                        O(k,1) = t; O(k,2) = r;
                    end
                end
                for t = 1:m
                    for r = 1:m
                        k = k+1;
                        EB(:,k) = orthopolyeval(C(t,:,j),X(1:N,j),t) .* ...
                            orthopolyeval(C(r,:,z),X(1:N,z),r);
                        O(k,1) = t; O(k,2) = r;
                    end
                end
                for i1 = 1:m
                    for i2 = 1:m
                        for i3 = 1:m
                            k = k+1;
                            EB(:,k) = orthopolyeval(C(i1,:,i),X(1:N,i),i1) .* ...
                                orthopolyeval(C(i2,:,j),X(1:N,j),i2) .* ...
                                orthopolyeval(C(i3,:,z),X(1:N,z),i3);
                            O(k,1) = i1; O(k,2) = i2; O(k,3) = i3;
                        end
                    end
                end
            end
        end
    end
end
return