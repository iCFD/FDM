%% %####################################################################%%%
%%%                   Print Mesh Refiment Test Results                  %%%
%%%#####################################################################%%%
clear all;
close all;           %remove items from workspace
clc;                 %clear command window
format compact       %clear empty rows

%% %####################################################################%%%
%%%                   Default Graphics parameters                       %%%
%%%#####################################################################%%%
% set(0,'defaultTextInterpreter','none') 
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',7.75)

%-------------------------------------------------------------------------
color=['k','r','g','c','m','b','y','w'];
lines={'-',':','--','-.','-.','none'};
mark=['s','+','o','x','v','none'];

%% %####################################################################%%% 
%%%               List of values for labeling                           %%%
%%%#####################################################################%%%

% norm error in properties
label_rho='$\rho/\rho_\infty$';
label_u='$ux/C_\infty$';
label_t='$\theta/\theta_\infty$';

% axis for L1 norm
label_y1='$\log\,\|e\|_1$';
label_x1='$\log\,\Delta{x}$';

% axis for L\inf norm
label_y2='$\log\,\|e\|_\infty$';
label_x2='$\log\,\Delta{x}$';

% for both norms
label_y3='$\log\,\|error\|_q$';
label_x3='$\log\,\Delta{x}$';

% Time intergration methods
%label_rk1='IMEX-SSP33';
label_rk1='IMEX3';
label_rk2='ARK3';
label_rk3='ARK4';

% Discretization methods
label_adv1='WENO3';
label_adv2='WENO5';
label_adv3='WENO7';

% norms name
label_Norm1 = '$\|e\|_1$';
label_Norm2 = '$\|e\|_2$';
label_Norm3 = '$\|e\|_\infty$';

% Result's folders names
label_result1='OOA-FDM-IMEX-SSP3';
label_result2='OOA-FDM-ARK3';
label_result3='OOA-FDM-ARK4';

%% %####################################################################%%%
%%%                    Save Figures to this path:                       %%%
%%%#####################################################################%%%
%path = '/home/manuel/Dropbox/Apps/Texpad/JCP2014/FDM_IMEX_AccuracyTest/';
path = '/home/manuel/Dropbox/Apps/Texpad/ECFD2014/FDM_IMEX_AccuracyTest/';

%% %####################################################################%%%
%%%                      Batch of Results to plot                       %%%
%%%#####################################################################%%%
results = {label_result1,label_result2,label_result3};
timeInt = {label_rk1,label_rk2,label_rk3};

%% %####################################################################%%%
%%%                 Defaul data from mesh refiment test                 %%%
%%%#####################################################################%%%

% number of nodes in every test:
nodes = [20 40 80 160 320]'; 

% cell size in every test:
dx = 1./(nodes-1);
       
%% %####################################################################%%%
%%%               Traditional Discrete Ordinate Method                  %%%
%%%#####################################################################%%%

for r = 1:3
    
folder = results{r};
load([folder,'/Norm.mat' ],'Norm' );
load([folder,'/OOA.mat'  ],'OOA'  );
load([folder,'/Stats.mat'],'Stats');

% find number of Methods and RKs schemes tested
rk = 1; %size(Norm,4); % one at the time!
mth = size(Norm,3);

%% %####################################################################%%%
%%%                 Print Norm L1 and L\infty results                   %%%
%%%#####################################################################%%%

ms=5;   % local marker size

% Set figure size and position
nameL1 = [results{r},'Lq_Norms'];
graph01=figure(r);
set(graph01, 'Position', [200 200 600 400]) 
set(gca,'YScale','log'); xlim([0.0001,0.1])
set(gca,'XScale','log'); ylim([1e-12,10])

% Create figure for Norm L1 for Density
hold on;
for irk = 1:rk
    for imth = 1:mth
        plot(dx,Norm(:,1,imth,1),[color(1),lines{1},mark(imth)],...
            'MarkerEdgeColor',color(1),...
            'MarkerSize',ms,...
            'MarkerFaceColor',color(8));
        plot(dx,Norm(:,2,imth,1),[color(2),lines{4},mark(imth)],...
            'MarkerEdgeColor',color(1),...
            'MarkerSize',ms,...
            'MarkerFaceColor',color(8));%,...
            %'LineWidth',1.5);
    end
end
% Draw reference triangles!
switch r % result#
    case {1,2}
        % Lower Triangle
        x1=1E-2; y1=1E-8; b1=2E-2; h1=3E-7;
        x_triag1 = [x1,x1+b1,x1+b1];
        y_triag1 = [y1,  y1 ,h1-y1];
        fill(x_triag1,y_triag1,'w','EdgeColor','r','LineWidth',1.5)
        % Upper Triangle
        x2=3E-3; y2=2E-2; b2=5E-3; h2=1E-1;
        x_triag2 = [x2,  x2 ,b2+x2];
        y_triag2 = [y2,y2+h2,y2+h2];
        fill(x_triag2,y_triag2,'w','EdgeColor','r','LineWidth',1.5)
        % Create textbox
        annotation(graph01,'textbox',...
            [0.50 0.75 0.023 0.046],...
            'String',{'3'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.56 0.79 0.023 0.046],...
            'String',{'1'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.71 0.31 0.023 0.046],...
            'String',{'1'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.78 0.38 0.023 0.046],...
            'String',{'3'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
    case {3}
        % Lower Triangle
        x1=1E-2; y1=1E-10; b1=2E-2; h1=1E-7;
        x_triag1 = [x1,x1+b1,x1+b1];
        y_triag1 = [y1,  y1 ,h1-y1];
        fill(x_triag1,y_triag1,'w','EdgeColor','r','LineWidth',1.5)
        % central Triangle
        x3=3E-3; y3=1E-7; b3=5E-3; h3=5E-6;
        x_triag3 = [x3,  x3 ,b3+x3];
        y_triag3 = [y3,y3+h3,y3+h3];
        fill(x_triag3,y_triag3,'w','EdgeColor','r','LineWidth',1.5)
        % Upper Triangle
        x2=3E-3; y2=2E-2; b2=5E-3; h2=1E-1;
        x_triag2 = [x2,  x2 ,b2+x2];
        y_triag2 = [y2,y2+h2,y2+h2];
        fill(x_triag2,y_triag2,'w','EdgeColor','r','LineWidth',1.5)
        annotation(graph01,'textbox',...
            [0.50 0.75 0.023 0.046],...
            'String',{'3'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.56 0.79 0.023 0.046],...
            'String',{'1'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.50 0.45 0.023 0.046],...
            'String',{'5'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.56 0.52 0.023 0.046],...
            'String',{'1'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.78 0.29 0.023 0.046],...
            'String',{'6'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
        % Create textbox
        annotation(graph01,'textbox',...
            [0.71 0.19 0.023 0.046],...
            'String',{'1'},...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'LineWidth',1);
end
hold off;

%grid on;

% Figure's Information
%title('Nodal DG, Norm error is $\rho(x)$, $\varepsilon = 1\times{10}^{-5}$');
xlabel(label_x3);
ylabel(label_y3);

% Figures details
leg_graph01=legend(...
    [label_Norm1,'-',label_adv1,'-',timeInt{r}],...
    [label_Norm3,'-',label_adv1,'-',timeInt{r}],...
    [label_Norm1,'-',label_adv2,'-',timeInt{r}],...
    [label_Norm3,'-',label_adv2,'-',timeInt{r}],...
    [label_Norm1,'-',label_adv3,'-',timeInt{r}],...
    [label_Norm3,'-',label_adv3,'-',timeInt{r}]);
set(leg_graph01,'Location','West','FontSize',14); legend('boxoff');

% Print Figure
print('-depsc',[path,nameL1,'.eps']);

end