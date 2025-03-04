function [w_upd,extB_upd,M_new] = HDMR_EXT_opt_order(w,extB,Y,f0, ...
    ho,M,maxorder,K_1,K_2,K_3,RMSEtol)
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
% This function computes optimum polynomial order based on RMSE           %
%                                                                         %
%  SYNOPSIS                                                               %
%   [w_upd,extB_upd,M_new] = HDMR_EXT_opt_order(w,extB,Y,f0, ...          %
%       ho,M,maxorder,K_1,K_2,K_3,RMSEtol)                                %
%                                                                         %
%  © Written by Jasper A. Vrugt & A. Sahin                                %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Make sure that w is a row vector
w = w(:); w = w';
% Size
[N,d] = size(extB);
% Allocation
Y_temp = zeros(N,d); Y_1 = zeros(N,M); Y_2 = Y_1; Y_3 = Y_1;
% Compute temporary matrix
for i = 1:N
    Y_temp(i,:) = extB(i,:) .* w;
end
% Built 1st Order Comp. Func. by ascending order
for i = 1:M
    idx1 = find(sum(ho(1:K_1,:),2) <= i);
    Y_1(:,i) = sum(Y_temp(:,idx1),2);
end
if maxorder > 1
    % Built 2nd Order Comp. Func. by ascending order
    for i = 1:M
        idx2_1 = sum(ho(K_1+1:K_1+K_2,1),2) <= i;
        idx2_2 = sum(ho(K_1+1:K_1+K_2,1),2) <= i;
        idx2 = find((idx2_1 + idx2_2) == 2) + K_1;
        Y_2(:,i) = sum(Y_temp(:,idx2),2);
    end
    if maxorder == 3
        % Built 3rd Order Comp. Func. by ascending order
        for i = 1:M
            idx3_1 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,1),2) <= i;
            idx3_2 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,2),2) <= i;
            idx3_3 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,3),2) <= i;
            idx3 = find((idx3_1 + idx3_2 + idx3_3) == 3) + K_1 + K_2;
            Y_3(:,i) = sum(Y_temp(:,idx3),2);
        end
    end
end
% Now, compute RMSE by order M
for i = 1:M
    if maxorder == 1
        Yt = Y_1(:,i);
    end
    if maxorder == 2
        Yt = Y_1(:,i) + Y_2(:,i);
    end
    if maxorder == 3
        Yt = Y_1(:,i) + Y_2(:,i) + Y_3(:,i);
    end
    RMSE(i) = sqrt(sum((Y - f0 - Yt).^2)/N);
end

% Determine optimum new order
idx = RMSE < RMSEtol;
if sum(idx) == 0
    M_new = find(RMSE == min(RMSE));
else
    idx = double(idx);
    idx(idx==0) = nan;
    M_new = min(find(idx == 1));
end
% Update extB and w
idx1 = find(sum(ho(1:K_1,:),2) <= M_new);
id_tot = idx1';
if maxorder > 1
    idx2_1 = sum(ho(K_1+1:K_1+K_2,1),2) <= M_new;
    idx2_2 = sum(ho(K_1+1:K_1+K_2,1),2) <= M_new;
    idx2 = find((idx2_1 + idx2_2) == 2) + K_1;
    id_tot = [id_tot, idx2'];
    if maxorder == 3
        idx3_1 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,1),2) <= M_new;
        idx3_2 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,2),2) <= M_new;
        idx3_3 = sum(ho(K_1+K_2+1:K_1+K_2+K_3,3),2) <= M_new;
        idx3 = find((idx3_1 + idx3_2 + idx3_3) == 3) + K_1 + K_2;
        id_tot = [id_tot, idx3'];
    end
end
extB_upd = zeros(N,d);
extB_upd(:,id_tot) = extB(:,id_tot);
w_upd = zeros(1,d);
w_upd(1,id_tot) = w(1,id_tot);


return