function qIR_MainFigure(LW,LETTERS,color,id0,Ny,T,days_per_year,t_begin_c, ...
      time,x,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,E0T,E1T,xmin,xmax, ...
      Temergence,RelativeGain)

      % MainFigure - Generates a 3x3 subplot figure with fitness functions and population data
      %
      % Inputs:
      %   LW          - Line width
      %   LETTERS     - Cell array of letters for subplot labels
      %   color       - Color matrix
      %   id0         - Index for initial condition
      %   Ny          - Number of years of insecticide usage
      %   T           - Final time
      %   t_begin_c   - Starting time of insecticide usage
      %   time        - Time vector
      %   x           - Phenotype vector
      %   Theta_i     - Fitness function of unexposed/exposed Adult Mosquitoes
      %   AA0, AA1    - AFM population matrices (a,x,t)
      %   E0, E1      - Egg population matrices (x,t)
      %   AA0T, AA1T  - Total unexposed/exposed AFM populations over time
      %   E0T, E1T    - Total egg populations over time
      %   Temergence  - Emergence time
      %   RelativeGain- Relative gain value


    %main figure
    figure
    set(gcf,'position',[100,100,1200,900]) %3eCoord=largeur, 4e=hauteur
    axes ('fontsize',15)

    %fitness functions
    ax=subplot(3,3,1.4);
    GrNum=1;
    yyaxis left 
    plot(x,Theta0,'LineWidth',LW,'color',color(1,:),'linestyle','-');
    xlabel('Phenotype ($x$)','Interpreter','latex','fontsize',15);
    xlim([xmin,xmax]);
    ylabel('$\Theta_0$','Interpreter','latex','fontsize',13);
    hold on
    yyaxis right 
    plot(x,Theta1,'LineWidth',LW,'color',color(2,:),'linestyle','--');
    ylabel('$\Theta_1$','Interpreter','latex','fontsize',13);
    hold off
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')} Fitness functions'], ...
        'interpreter','latex');
     ax.YAxis(1).Color = color(1,:);
     ax.YAxis(2).Color = color(2,:);


    %Rgain and Temg
    GrNum=GrNum+1;
    subplot(3,3,2.6);
    y=AA0T(id0)+AA1T(id0);
    z=AA0T+AA1T;
    plot([0 T],[y y],'LineWidth',LW,'color',color(1,:),'linestyle','-');
    xlim([0 T]);
    ylim([min(z) max(z)]);
    xlabel('Time $t$ (year)','Interpreter','latex','fontsize',15);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('AFMs ($A_0+A_1$)','Interpreter','latex','fontsize',13);
    hold on
    plot(time,z,'LineWidth',LW,'color',color(2,:),'linestyle','--');
    TemergPlot=Temergence+t_begin_c;
    plot([TemergPlot TemergPlot],[min(z) max(z)],'LineWidth',.5,'color','k','linestyle',':');
    hold off
    title(['\fontsize{13}{0}\selectfont' '\quad' '\textbf{(' LETTERS(GrNum) ')}' ...
           '\fontsize{13}{0}\selectfont$\:{\rm r}_{\rm gain}=$' num2str(RelativeGain,2) ...
           '\ ; \fontsize{13}{0}\selectfont$\:{\rm T}_{\rm emg}=$' num2str(Temergence/days_per_year,3)], ...
           'interpreter','latex');
    legend('$c=0$','$c>0$','Interpreter','latex','location','best',...
        'Orientation','vertical','fontsize',13);
    legend boxoff 


    %AFM population
    GrNum=GrNum+1;
    subplot(3,3,4);
    surf(time,x,AA0','FaceColor','interp', 'EdgeColor','none', ...
          'FaceLighting','gouraud');
    axis tight
    view(-136,60);
    camlight left
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13,'rotation',-30);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',13,'rotation',35);
    zlabel('$A_0(t,x)$','Interpreter','latex','fontsize',15);
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');

    GrNum=GrNum+1;
    subplot(3,3,5);
    surf(time,x,AA1','FaceColor','interp', 'EdgeColor','none', ...
          'FaceLighting','gouraud');
    axis tight
    view(-136,60);
    camlight left
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13,'rotation',-30);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',13,'rotation',35);
    zlabel('$A_1(t,x)$','Interpreter','latex','fontsize',15);
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');

    GrNum=GrNum+1;
    ax=subplot(3,3,6);
    yyaxis left 
    plot(time,AA0T,'LineWidth',LW,'color',color(1,:),'linestyle','-');
    xlim([0 T]);
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('$A_0(t)$','Interpreter','latex','fontsize',15);
    hold on
    yyaxis right 
    plot(time,AA1T,'LineWidth',LW,'color',color(2,:),'linestyle','--');
    ylabel('$A_1(t)$','Interpreter','latex','fontsize',15);
    hold off
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');
    ax.YAxis(1).Color = color(1,:);
    ax.YAxis(2).Color = color(2,:);

    
    %eggs population
    GrNum=GrNum+1;
    subplot(3,3,7);
    surf(time,x,E0','FaceColor','interp', 'EdgeColor','none', ...
          'FaceLighting','gouraud');
    axis tight
    view(-136,60);
    camlight left
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13,'rotation',-30);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',13,'rotation',35);
    zlabel('$E_0(t,x)$','Interpreter','latex','fontsize',15);
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');

    GrNum=GrNum+1;
    subplot(3,3,8);
    surf(time,x,E1','FaceColor','interp', 'EdgeColor','none', ...
          'FaceLighting','gouraud');
    axis tight
    view(-136,60);
    camlight left
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13,'rotation',-30);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',13,'rotation',35);
    zlabel('$E_1(t,x)$','Interpreter','latex','fontsize',15);
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');

    GrNum=GrNum+1;
    ax=subplot(3,3,9);
    yyaxis left 
    plot(time,E0T,'LineWidth',LW,'color',color(1,:),'linestyle','-')
    xlim([0 T])
    xlabel('Time ($t$)','Interpreter','latex','fontsize',13);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('$E_0(t)$','Interpreter','latex','fontsize',15);
    hold on
    yyaxis right 
    plot(time,E1T,'LineWidth',LW,'color',color(2,:),'linestyle','--');
    ylabel('$E_1(t)$','Interpreter','latex','fontsize',15);
    hold off
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
        'interpreter','latex');
    ax.YAxis(1).Color = color(1,:);
    ax.YAxis(2).Color = color(2,:);

end
