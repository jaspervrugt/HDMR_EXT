function [SA,SA_sig,Fx] = HDMR_EXT_end(SA,nterms,p0,RMSE,select,K, ...
    C2,C3,n1,n2,n3,n,n_ns)
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
% Prepares the return arguments of HDMR_EXT                               %
%                                                                         %
%  SYNOPSIS                                                               %
%   SA,SA_sig,Fx] = HDMR_EXT_end(SA,nterms,p0,RMSE,select,K, ...          %
%       C2,C3,n1,n2,n3,n,n_ns)                                            %
%                                                                         %
%  REFERENCE                                                              %
%   Genyuan Li et al., Journal of Physical Chemistry A., V. 114 (19),     %
%       pp. 6022-6032, 2010.                                              %
%                                                                         %
%  � Written by Jasper A. Vrugt, Sept. 9, 2017                            %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Skip line with HDMR print statement about progress
fprintf('\n');
% Set alfa [significance level] for boostrap confidence intervals
alfa = 0.05;

% ----------------------------------------------------------------------- %
%          Now compile table, SI, with all sensitivity results           |%
% ----------------------------------------------------------------------- %

% Take sum of all selected terms of all K bootstrap trials
f_ret = sum(select,2);
% Now construct two tables with results
[SA,Fx] = construct_tables(SA,nterms,p0,RMSE,alfa,f_ret,K,C2,C3,n1, ...
    n2,n3,n,n_ns);
% Compute some variables - used to write tables
nr = size(SA,1); idx_sig = [ 1 ; 1 + find( f_ret > 0) ; nr ]; 
SA_sig = SA(idx_sig,1:13);

% ----------------------------------------------------------------------- %
%                   Print a Table with emulator results                   %
% ----------------------------------------------------------------------- %

% Open file for readme.txt
fid = fopen('HDMR_EXT_results.txt','w+');
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
fprintf(fid,'|   RRRRR     EEEEEEEE   SSSSSSS  UU    UU  LL        TTTTTTTT  SSSSSSSS  |\n');
fprintf(fid,'|   RRRRRR    EEEEEEEE  SSSSSSS   UU    UU  LL        TTTTTTTT  SSSSSSS   |\n');
fprintf(fid,'|   RR   RR   EE        SS        UU    UU  LL        TT TT TT  SS        |\n');
fprintf(fid,'|   RR   RR   EEEEE     SSSSSS    UU    UU  LL           TT     SSSSSS    |\n');
fprintf(fid,'|   RRRRRR    EEEEE     SSSSSSS   UU    UU  LL           TT     SSSSSSS   |\n');
fprintf(fid,'|   RR  RR    EE              SS  UUU  UUU  LL    LL     TT           SS  |\n');
fprintf(fid,'|   RR   RR   EEEEEEEE   SSSSSSS  UUUUUUUU  LLLLLLLL     TT      SSSSSSS  |\n');
fprintf(fid,'|   RR    RR  EEEEEEEE  SSSSSSS    UUUUUU   LLLLLLLL     TT     SSSSSSS   |\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|-------------------------------------------------------------------------|\n');
fprintf(fid,'\n');
% Now print table with content of cell matrix B
if K == 1
    fprintf(fid,['Table 1: Properties of HDMR emulator, y = f(x), ' ...
        'for training data set (no bootstrapping)\n']);
elseif K > 1
    fprintf(fid,['Table 1: Properties of HDMR emulator, y = f(x), ' ...
        'for randomized training data set of %d bootstrap trials\n'],K);
end
fprintf(fid,'===============================================\n');
fprintf(fid,'Emulator    # terms   # coefs.   RMSE     ns   \n');
fprintf(fid,'-----------------------------------------------\n');
% Define format of printing for table
fmt_1 = ('   %-7d\t %3d        %3d      %5.3f    %3d\n');
for k = 2:K+1
    % print according to prespecified format
    fprintf(fid,fmt_1,str2double(Fx{k,1}),Fx{k,2},Fx{k,3}, ...
        str2double(Fx{k,4}),Fx{k,5});
end
fprintf(fid,'===============================================\n');
fprintf(fid,'\n');

% ----------------------------------------------------------------------- %
%                  Print a table with sensitivity results                 %
% ----------------------------------------------------------------------- %

% Print to file the results of the HDMR analysis
for pr_tab = 1:2
    % Note: Statement \261 prints "�" symbol
    fprintf(fid,'\n'); fprintf(fid,'\n');
    if ( pr_tab == 1 )      % --> print full table with all terms included
        idx_table = 1:nr; table = SA;
    elseif ( pr_tab == 2 )  % --> print only significant terms to table
        idx_table = idx_sig; nr = numel(idx_sig); table = SA_sig;
    end
    % Print table header in HDMR_results.txt
    if K == 1
        if pr_tab == 1
            fprintf(fid,['Table 2: HDMR_EXT results for all model ' ...
                'components using all X and Y data (no bootstrapping)\n']);
        else
            fprintf(fid,['Table 3: HDMR_EXT results for significant model ' ...
                'components only using all X and Y data ' ...
                '(no bootstrapping)\n']);
        end
    elseif K > 1
        if pr_tab == 1
            fprintf(fid,['Table 2: HDMR_EXT results for all model ' ...
                'components using %d bootstrap trials\n'],K);
        else
            fprintf(fid,['Table 3: HDMR_EXT results for significant model ' ...
                'components only using %d bootstrap trials\n'],K);
        end
    end
    % Now print table
    fprintf(fid,['====================================================' ...
        '==========================================\n']);
    fprintf(fid,['                                               ( = J' ...
        'OINT DETERMINATION )                      \n']);
    fprintf(fid,['                  ----------------------------------' ...
        '------------------------------------------\n']);
    fprintf(fid,['Term       Order      S^a           S^b           S ' ...
        '            ST         V(i)/V(Y)   #select\n']);
    fprintf(fid,['----------------------------------------------------' ...
        '------------------------------------------\n']);
    % Define format of printing for table
    if ( K == 1 )
        fmt_1 = (['%-11s  %1d \t %5.2f (-----) %5.2f (-----) %5.2f ' ...
            '(-----) %5.2f (-----) %5.2f (-----)   %-3d\n']);
        fmt_2 = (['%-11s  %1d \t %5.2f (-----) %5.2f (-----) %5.2f ' ...
            '(-----)       (-----) %5.2f (-----)   %-3d\n']);
        for r = 2:nr-1
            if idx_table(r-1) <= n1
                % isolate data
                data = str2double(table(r,[3 5 7 9 11]))';
                % print according to prespecified format
                fprintf(fid,fmt_1,char(table{r,1}),table{r,2}, ...
                    data(1:5),table{r,13});
            else
                % isolate data
                data = str2double(table(r,[3 5 7 11]))';
                % print according to prespecified format
                fprintf(fid,fmt_2,char(table{r,1}),table{r,2}, ...
                    data(1:4),table{r,13});
            end
        end
    elseif ( K > 1 )
        fmt_1 = (['%-11s  %1d \t %5.2f (\261%.2f) %5.2f (\261%.2f) ' ...
            '%5.2f (\261%.2f) %5.2f (\261%.2f) %5.2f (\261%.2f)   %-3d\n']);
        fmt_2 = (['%-11s  %1d \t %5.2f (\261%.2f) %5.2f (\261%.2f) ' ...
            '%5.2f (\261%.2f)               %5.2f (\261%.2f)   %-3d\n']);
        for r = 2:nr-1
            if idx_table(r-1) <= n1
                % isolate data
                data = str2double(table(r,3:12))';
                % print according to prespecified format
                fprintf(fid,fmt_1,char(table{r,1}),table{r,2}, ...
                    data(1:10),table{r,13});
            else
                % isolate data
                data = str2double(table(r,[3:8 11:12]))';
                % print according to prespecified format
                fprintf(fid,fmt_2,char(table{r,1}),table{r,2}, ...
                    data(1:8),table{r,13});
            end
        end
    end
    fprintf(fid,['----------------------------------------------------' ...
        '------------------------------------------\n']);
    if ( K == 1 )
        %5.2f
        fmt = (['%-11s  %1d \t %5.2f (-----) %5.2f (-----) %5.2f ' ...
            '(-----)       (-----) %5.2f (-----)\n']);
        % isolate data
        data = str2double(table(nr,[3 5 7 11]))'; fprintf(fid,fmt, ...
            char(table{nr,1}),[],data(1:4));
    else
        %5.2f (-----)
        fmt = (['%-11s  %1d \t %5.2f (\261%.2f) %5.2f (\261%.2f) ' ...
            '%5.2f (\261%.2f)               %5.2f (\261%.2f) \n']);
        % isolate data
        data = str2double(table(nr,[3:8 11:12]))'; fprintf(fid,fmt, ...
            char(table{nr,1}),[],data(1:8));
    end
    fprintf(fid,['====================================================' ...
        '==========================================\n']);
    fprintf(fid,[' S^a: Structural sensitivity index of individual ' ...
        'terms\n']);
    fprintf(fid,[' S^b: Correlative sensitivity index of individual ' ...
        'terms\n']);
    fprintf(fid,' S: Total sensitivity index of individual terms\n');
    fprintf(fid,' ST: Total sensitivity index\n');
    if K == 1
        fprintf(fid,[' (--): Cannot compute confidence intervals of ' ...
            'listed statistics with K = 1\n']);
    elseif K > 1
        fprintf(fid,[' (\261 ): %2d%% confidence intervals derived ' ...
            'from bootstrapping\n'],100*(1-alfa));
    end
    fprintf(fid,[' V(i)/V_Y: Relative contribution of each term to ' ...
        'model output variance ( = var(Y) ) \n']);
    if K == 1
        fprintf(fid,[' #select: 0 (if term is insignificant) and 1 ' ...
            '(significant term)\n']);
    elseif K > 1
        fprintf(fid,[' #select: Number of bootstrap trials that ' ...
            'identifies respective term as significant\n']);
    end
    fprintf(fid,'\n');
end
% close the RMSE_output.txt file
fclose('all');
% Now open in text files in editor (not on linux)
if ispc || ismac
    edit HDMR_EXT_settings.txt; edit HDMR_EXT_results.txt;
end