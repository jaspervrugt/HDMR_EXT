function HDMR_EXT_plot(SA_sig,Fx,y,Y_e,select,p0,n_ns,id,R,K)
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
% Plots the results of HDMR_EXT                                           %
%                                                                         %
%  SYNOPSIS                                                               %
%   HDMR_EXT_plot(SA_sig,Fx,y,Y_e,select,p0,n_ns,id,R,K)                  %
%                                                                         %
%  © Written by Jasper A. Vrugt, Sept. 9, 2017                            %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Print wait statement to the screen
fprintf('HDMR_EXT PLOTTING: PLEASE WAIT ...');
% Define name of program
n_program = 'HDMR_EXT'; 
% Define name of figures file
file_name = [n_program,'_figures.pdf'];
% Define names of tables
table1_name = ['Table_1_',n_program,'.pdf']; 
table2_name = ['Table_2_',n_program,'.pdf']; 
% Determine the entries without pdf
id_pdf = 1:(strfind(table1_name,'pdf') - 2);

% Determine screen size
scr_z = get(0,'ScreenSize'); 
% Multiplier, x and y: axis
x_mult = scr_z(3)/1920; y_mult = scr_z(4)/1080;
% Multiplier, text
t_mult = min(x_mult,y_mult);
% Define fontsize for figures
fontsize_xylabel = 18 * t_mult;
fontsize_axis = 16 * t_mult;
fontsize_legend = 14 * t_mult;
fontsize_text = 14 * t_mult;
fontsize_table = 16 * t_mult;
fontsize_titlepage = 30 * t_mult;

% ----------------------------------------------------------------------- %
%                 Now plot empty figure for PDF file                      %
% ----------------------------------------------------------------------- %

figure('units','normalized','outerposition',[0 0 1 1],'name','first page');
plot([],[],'ro'); axis([0 1 0 1]); set(gcf,'color','w'); 
set(gca,'XColor','w','YColor','w');
%title('Visual results of HDMR toolbox','fontsize',20,'interpreter','latex');
text(0.3*x_mult,0.6*y_mult,'Visual results of HDMR$_{\rm EXT}$ toolbox', ...
    'fontsize',fontsize_titlepage,'interpreter','latex');
% Now about Tables
text(0.3*x_mult,0.5*y_mult,'$\;\;\;\;\;$GUI Tables may not print well ', ...
    'fontsize', fontsize_titlepage,'interpreter','latex');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                 Now plot results of emulator of K trials                %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% First plot scatter plot with results of each emulator
Y_plot = linspace(min(y),max(y));
% plot emulators in subplots with r rows and c columns
row = 2; col = 5;
% Now loop over each emulator and make scatter plot
for k = 1 : K
    % If true: open new figure and reset nr of columns to one
    if (mod(k,row*col) == 1) || (k == 1)
        figure('units','normalized','outerposition',[0 0 1 1],...
            'name',['HDMR_EXT: Bivariate scatter plots actual and ' ...
            'emulated y values of training data set']); c = 1;
    end
    % determine appropriate index of figure
    if c <= col
        r = 1; plot_c = c;
    else
        r = 2; plot_c = c - col;
    end
    % select right Y indices
    id_R = id(1:R,k);
    % Define the location of the figure
    ax1 = axes('units','inches');
    % define new axis position
    axpos1 = [0.9 + (plot_c - 1) * 3.75 6 - (r-1) * 4.75 3 3]; 
    % scale in x and y
    axpos1([1 3]) = axpos1([1 3]) * x_mult;
    axpos1([2 4]) = axpos1([2 4]) * y_mult;
    % End scale in x and y
    set(ax1,'position',axpos1);
    % Now plot the training data
    plot(y(id_R),Y_e(1:R,k),'r+','markersize',12,'linewidth',2); 
    hold on; axis square; axis tight;
    % Now plot the 1:1 Line in gray
    plot(Y_plot,Y_plot,'k-','color',[0.5 0.5 0.5],'linewidth',2); 
    axis tight;
    % Add labels
    xlabel('$y$','fontsize',fontsize_xylabel,'interpreter','latex');
    if mod(k,col) == 1
        ylabel('$y = f(\textbf{x})$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
    end
    set(gca,'fontsize',fontsize_axis);
    % add legend only in top-left plot
    if ( c == 1 )
        leg = legend({'{\color{red} Training}',...
            '{\color{gray} 1:1 Line}'},'location', ...
            'southeast'); %,'interpreter','latex');
        set(leg,'linewidth',3); set(leg,'fontsize',fontsize_legend); 
        legend boxoff;
    end
    % Calculate training statistics
    RMSE_train = sqrt(1/R * sum((y(id_R) - Y_e(1:R,k)).^2));
    RMSE_1 = num2str(RMSE_train); RMSE_1 = RMSE_1(1:min(numel(RMSE_1),5));
    title(strcat('Bootstrap trial:',{' '},num2str(k)), ...
        'fontweight','bold','interpreter','latex');
    % Now print results/statistics
    a = axis;
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.92*(a(4)-a(3)), ...
        strcat('$\# terms.:',{' '},...
        num2str(sum(select(:,k))),'$'),'interpreter','latex', ...
        'fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.85*(a(4)-a(3)), ...
        strcat('$\# coef.:',{' '},...
        num2str(p0(k)),'$'),'interpreter','latex', ...
        'fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.78*(a(4)-a(3)), ...
        strcat('$\# nonsgl.:',{' '},...
        num2str(n_ns(k)),'$'),'interpreter','latex', ...
        'fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.70*(a(4)-a(3)), ...
        strcat('${\rm RMSE}_{\rm T}:',{' '},RMSE_1,'$'),...
        'interpreter','latex','fontsize',fontsize_text);
    % update column number
    c = c + 1;
end
% Get figure handles
figHandles = flipud(findall(0,'Type','figure'));
for zz = 1:numel(figHandles)
    figure(figHandles(zz)); set(gcf,'color','w');
    switch zz
        case 1
            exportgraphics(figHandles(1),file_name);
        otherwise
            % Append to existing figure
            exportgraphics(figHandles(zz),file_name,'Append',true);
    end
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                     Now plot Fx to table in figure                      %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

col_names = Fx(1,1:5); n_col = numel(col_names);
for j = 1:n_col
    col_names{j} = sprintf(['<html><center /><font size=5> ' ...
        '%s </font></html>'],char(col_names(j)));
end
row_names = Fx(2:K+1,1);
for i = 1:K
    row_names{i} = sprintf(['<html><font size=5> ' ...
        '%s </font></html>'],char(row_names{i,1}));
end
% now justify content of each cell
for i = 2:K+1
    j = 1; Fx(i,j) = strcat({'         '},num2str(Fx{i,j}),{' '});  
    for j = 2:3
        Fx(i,j) = strcat({'       '},num2str(Fx{i,j}),{' '}); 
    end
    j = 4; Fx(i,j) = strcat({'     '},num2str(Fx{i,j}),{' '});
    j = 5; Fx(i,j) = strcat({'     '},num2str(Fx{i,j}),{' '});
end
% approximate height of table
ht = min(0.8,(K - 1) * 0.04 + 0.20);
% open new figure called Table 1: Emulator performance training data set
axpos_fig = [ 0.15 0.05 0.48 ht ]; 
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;

% open new figure called Table 1: Emulator performance training data set
fig = figure('units','normalized','name',['Table 1: Emulator performance ' ...
    'training data set'],'numbertitle','off',...
    'outerposition',axpos_fig);
% plot table
axpos_table = [ 0.1 0.04 0.82 0.6 ];
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y
tbl = uitable(fig,'units','normalized','position',axpos_table,...
   'Data',Fx(2:K+1,1:5),'ColumnName',col_names,'rowname', ...
   row_names,'fontsize',fontsize_table);
set(tbl,'columnWidth', {120})
% set row_header/cell_size
set_row_cells(tbl,col_names,row_names);

% Add title
axpos_fig = [ 0.07 0.7 (1 - 2*0.07) 0.25 ];
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;
% end scale in x and y
if K == 1
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_fig , 'string', ...
        ['Table 1: Properties of HDMR_EXT emulator, y = f(x), ' ...
        'for training data set (no bootstrapping)'],...
        'fontsize',fontsize_xylabel);
elseif K > 1
    str = strcat(['Table 1: Properties of HDMR_EXT emulator, ' ...
        'y = f(x), for randomized training data set of'],...
        {' '},num2str(K),{' '},'bootstrap trials');
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_fig, 'string', str, 'fontsize',fontsize_xylabel);
end
% Save as a PDF
print(gcf,table1_name(id_pdf),'-dpdf','-fillpage'); 

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   Now plot SA_sig to table in figure                    %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Now print figure with Tabulated results ( similar to HDMR_results.txt )
n_col = size(SA_sig,2); col_names = SA_sig(1,2:n_col); 
for j = 1:n_col-1
    col_names{j} = sprintf(['<html><center /><font size=5> ' ...
        '%s </font></html>'],char(col_names(j)));
end
n_row = size(SA_sig,1); row_names = cell(n_row-1,1);
for i = 2:n_row
    row_names{i-1} = sprintf(['<html><font size=5> ' ...
        '%s </font></html>'],char(SA_sig{i,1}));
end
% now justify content of each cell
for i = 2:n_row
    j = 2;  SA_sig(i,j) = strcat({'        '},num2str(SA_sig{i,j}),{' '});
    for j = 3:12
        SA_sig(i,j) = strcat({'     '},SA_sig(i,j),{' '}); 
    end
    j = 13;  SA_sig(i,j) = strcat({'       '},num2str(SA_sig{i,j}),{' '});
end
% approximate height of table
ht = min(0.8, (n_row - 1) * 0.038 + 0.15);
% open new figure called Table 1: Emulator performance training data set
axpos_fig = [ 0.02 0.2 0.85 ht ]; 
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;
% end scale in x and y

% "Table 2: Variance-based decomposition and sensitivity coefficients"
figure('units','normalized','name',['Table 2: Variance-based ' ...
    'decomposition and sensitivity coefficients'],'numbertitle','off',...
    'outerposition',axpos_fig);
% plot table
axpos_table = [ 0.04 0.1 0.90 0.73 ];
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y
tbl = uitable('units','normalized','position',axpos_table,...
   'data',SA_sig(2:n_row,2:n_col),'columnname',col_names, ...
   'rowname',row_names,'fontsize',fontsize_table);
set(tbl,'columnWidth', {100,120,100,120,100,120,100,120, ...
    100,135,100,120,135,100})
% set row_header/cell_size
set_row_cells(tbl,col_names,row_names)

% Now print title
% get(gcf, 'Position');
axpos_table = [ 0.07 0.85 (1 - 2*0.07) 0.12 ];
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y
if K == 1
    uicontrol('units','normalized','style', 'text', ...
        'position', axpos_table , 'string', ...
        ['Table 2: HDMR_EXT results for significant model ' ...
        'components only using all X and Y data (no bootstrapping)'],...
        'fontsize',fontsize_xylabel);
elseif K > 1
    str = strcat(['Table 2: HDMR_EXT results for significant model ' ...
        'components only using'],...
        {' '},num2str(K),{' '},'bootstrap trials');
    uicontrol('units','normalized','style', 'text', 'position', ...
        axpos_table, 'string', str, 'fontsize',fontsize_xylabel);
end
% Save as a PDF
print(gcf,table2_name(id_pdf),'-dpdf','-fillpage'); 

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%               Save Figures & Tables to a single PDF file                %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Merge figures & tables
input_files = {file_name,table1_name,table2_name};
% Current directory
HDMR_dir = pwd; 
% Forward slash or backward slash?
if contains(HDMR_dir,'/'), slh = '/'; else, slh = '\'; end
% How many PDFs?
n_pdfs = numel(input_files);
% Cell string
input_pdfs = cell(1,n_pdfs);
% Full names of input files
for i = 1:n_pdfs
    input_pdfs(i) = {[HDMR_dir,slh,input_files{i}]};
end
% Name of output file
output_pdf = [HDMR_dir,slh,file_name]; 
% Benjamin Großmann (2021). Merge PDF-Documents 
% (https://www.mathworks.com/matlabcentral/fileexchange/89127- ...
%     merge-pdf-documents), MATLAB Central File Exchange. August 24, 2021.
memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
cellfun(@(f) merger.addSource(f), input_pdfs)
merger.setDestinationFileName(output_pdf)
merger.mergeDocuments(memSet)

% Delete Tables
delete(table1_name); delete(table2_name); 
% Open PDF document
open(file_name);

% Print wait statement to the screen
fprintf(' DONE\n');

end
