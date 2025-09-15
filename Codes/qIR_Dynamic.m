clc
close all

fpath = '/Users/kezetasylvere/Documents/Codes IRSpread';   %%%%%%%%%%%%%%%%

LW= 1.5;     %Line width

LETTERS='A':'Z';  %to number the panels in plots with many subplots

color= [0.255 0.411 0.882;...   %royal blue
        0.635 0.078 0.184;...   %marron 
        0.133 0.698 0.133;...   %forest green
        0.890 0.260 0.200];...  %vermillon (rouge)

%% SECTION 1: Set of fixed parameters

%time
Ny= 10;                          %number of years of insecticide usage
days_per_year=300;               %number of days per year
t_begin_c= days_per_year;        %starting time of insecticide usage
T= t_begin_c+Ny*days_per_year;   %Total time
dt= 1;                           %time step
time= 0:dt:T;                    %discretisation of time
Nt= length(time);                %number of time steps

%age
A= 30;                      %mosquito age limit
AvrMosquiLifeSpan= 21;      %average mosquito life span in days
MinMosquiAgeLayEggs= 3;     %minimum age for AFM to lay eggs in days
da= .1;                     %age step
age= 0:da:A;                %discretisation of age
Na= length(age);            %number of age steps

%resistance   
xmin= -0.5;         %minimum value of insecticide resistance rate
xmax= 1.5;          %maximum value of insecticide resistance rate
x0= 0;              %Reference sensitive  insecticide resistance level                          
x1= 1;              %Reference resistant insecticide resistance level    
VarJ0= 0.01;        %Variance de J0
VarJ1= 2*VarJ0;     %Variance de J1
dx= VarJ0/5;        %insecticide rate step (dx<VarJ0)
x= xmin:dx:xmax;    %discretisation of resistance level
Nx= length(x);      %number of resistance steps
    
%death rates parameters
ProbaSurviFullySensi= 10^-10;       %daily probability of surviving insecticide exposure
dd0= -log(ProbaSurviFullySensi);    %death rate of ref sensi strain due to insecticide
dd1= -log(1-ProbaSurviFullySensi);  %death rate of ref resis strain due to insecticide

%AFM laying eggs rate parameters
rrm= 200;          %maximum number of eggs laid by the ref sensi strain
rr0= 0.875*rrm;    %average number of eggs laid by the ref sensi strain
rr1= 0.75*rrm;     %average number of eggs laid by the ref resis strain
             
%natural death rate of eggs laid by unexposed/exposed AFM
mu0= 0.2*ones(1,Nx);             
mu1= mu0;                           

%hatching rate of eggs laid by unexposed/exposed AFM
gamma0= 0.75*ones(1,Nx);          
gamma1= gamma0;                      

%proportion of hatched eggs laid by AFM that reach adulthood
tau= 0.67*ones(1,Nx);

%other fixed parameters
EmergenceThreshold=10;   %in percentage
k=0.5;                   %exposant of the logistic function growth of eggs   

%% SECTION 2: Table generated to quantify the effect of C on Rgain and Temg
if (0)
    format long
    tic

    %parameters setting
    C_est_cte=1; %if 1 then exposure rate is constant
    DurCycle=50; %cycle duration
    T0=200; %time between two cycles
    cVal= zeros(1,10); %vector of exposure rates

    %discretization of exposure rates
    VectC=0.05:0.05:0.95; lc=length(VectC);

    %storage of inputs
    ColC = VectC';  % Pre-populate known values
    ColRgain=zeros(lc,1);
    ColTemg=zeros(lc,1);

    %begin experience plan
    plan_exp = [VectC', -111*ones(lc,2)];
    % row_number=0;
    % for cm= 0.05:0.05:0.95
    %     row_number=row_number+1;
    %     plan_exp(row_number,1)=cm;
    % end
    
    % ===== PROGRESS MONITORING SETUP =====
    fprintf('===== SIMULATION STARTING =====\n');
    fprintf('Total number of simulations to run: %d\n', lc);
    fprintf('Starting simulations...\n\n');
    
    %begin computation on experience plan
    for row_number=1:lc
        %display progress: current run / total runs
        fprintf('==Processing simulation %d/%d==\n', row_number, lc);
        
        %start timing for individual simulation
        tic_sim = tic;
        
        param=plan_exp(row_number,:);
        cm=param(1);
        [~,age,x,~,~,~,~,~,~,~,~,~,~,~,~,RelativeGain,Temergence]= ...
              qIR_MainFunction(C_est_cte,cm,VarJ0,VarJ1,time,T,t_begin_c, ...
              tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0,gamma1,mu0,mu1,rrm,rr0, ...
              rr1,dd0,dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs,cVal,k, ...
              EmergenceThreshold);

        ColRgain(row_number)=RelativeGain;
        ColTemg(row_number)=Temergence;
        
        %show completion of current simulation with time
        sim_time = toc(tic_sim);
        fprintf(['---> Simulation completed in %.2f seconds --- ' ...
              ' %.1f%% completed\n\n'], sim_time, (row_number/lc)*100);
    end
    
    % Multiple beeps to signal completion
    for i = 1:5
        beep
        pause(0.5)  % Half-second pause between beeps
    end

    %final completion message
    fprintf('===== ALL SIMULATIONS COMPLETED =====\n');
    fprintf('Total simulations run: %d\n', lc);
    
    %outputs table name
    tabName='TabRgainTemg.xlsx';

    %here we generate the contents of the simulation table
    Tab=table(ColC,ColRgain,ColTemg);
    writetable(Tab,fullfile(fpath,tabName));
    
    fprintf('Results saved to: %s\n', tabName);
    
    total_time = toc;
    fprintf('Total execution time: %.2f seconds (%.2f minutes)\n', total_time, total_time/60);
end

%% SECTION 3: Figure effect of C on Rgain and Temg
if (0)

    Data=readtable('/Users/kezetasylvere/Documents/Codes IRSpread/TabRgainTemg.xlsx');

    %vector of exposure rates
    VectC=0.05:0.05:0.95;                  lc=length(VectC);

    %extract Rgain and Temg
    yRgain=zeros(lc,1); 
    yTemg=zeros(lc,1);
    for j=1:lc
        yRgain(j)=Data.ColRgain(j);
        yTemg(j)=Data.ColTemg(j);
    end

    %draw figure
    figure;
    ax=gca;
    yyaxis left 
    plot(VectC,yRgain,'LineWidth',LW,'color',color(1,:),'linestyle','-')
    xlabel('Exposure rate ($c$)','Interpreter','latex','fontsize',18)
    ylabel('${\rm r}_{\rm gain}$','Interpreter','latex','fontsize',18)
    hold on
    yyaxis right 
    plot(VectC,yTemg,'LineWidth',LW,'color',color(2,:),'linestyle','--')
    ylabel('${\rm T}_{\rm emg}$','Interpreter','latex','fontsize',18)
    hold off

    ax.YAxis(1).Color = color(1,:);
    ax.YAxis(2).Color = color(2,:);

    %save the figure
    fileName1 = 'FigRgainTemg.pdf';
    exportgraphics(gcf, fullfile(fpath, fileName1), ...
    'ContentType', 'image', 'BackgroundColor','white');
end 

%% SECTION 4: Fitness functions, survival probabilities, death and egg-laying rates=fct(C_est_cte,cm,cVal,ProbaSurviFullySensi)
if (0)
    %setting parameters values
    C_est_cte=1;
    cm=0.5;
    cVal= zeros(1,Ny);
    ProbaSurviFullySensi=10^-10;
    dd0= -log(ProbaSurviFullySensi);   
    dd1= -log(1-ProbaSurviFullySensi);  

    %insecticide exposure rate
    id0=1+floor(t_begin_c/dt);    
    c= zeros(1,Nt);
    if C_est_cte
        for t=(id0+1):Nt
             c(t)=cm;
        end
    else
        for j=1:10
            Tdebut=t_begin_c+(j-1)*300;
            IdDebut=1+floor(Tdebut/dt);
            Tfin=Tdebut+300;
            IdFin=1+floor(Tfin/dt);
            for ll=IdDebut:IdFin
                c(ll)=cVal(j);
            end
        end
    end
    
    %run qIR_FunctionTheta
    [Theta0,Theta1,r0,r1,r0X,r1X,SurvProb0,SurvProb1,DeathInsecticide,d0,d1]=...
      qIR_FuncTheta(tau,age,x,Na,Nx,gamma0,gamma1,mu0,mu1,rrm,rr0,rr1,dd0,dd1, ...
      AvrMosquiLifeSpan,MinMosquiAgeLayEggs);

    %fitness functions
    figure
    yyaxis left 
    plot(x,Theta0,'LineWidth',LW,'color',color(1,:),'linestyle','-');
    xlabel('Phenotype ($x$)','Interpreter','latex','fontsize',15);
    xlim([-0.5 1.5]);
    ylabel('$\Theta_0$','Interpreter','latex','fontsize',15);
    hold on
    yyaxis right 
    plot(x,Theta1,'LineWidth',LW,'color',color(2,:),'linestyle','--');
    ylabel('$\Theta_1$','Interpreter','latex','fontsize',15);
    hold off
    title('Fitness functions','interpreter','latex','fontsize',15);
     ax.YAxis(1).Color = color(1,:);
     ax.YAxis(2).Color = color(2,:);
    %save the figure
    fileName2= 'Fitness functions.pdf';
    exportgraphics(gcf, fullfile(fpath, fileName2), ...
    'ContentType', 'image', 'BackgroundColor','white');


    %Exposure rate
    figure
    plot(time,c,'LineWidth',LW,'Color',color(4,:),'LineStyle','-.');
    xlabel('Time ($t$)','Interpreter','latex','fontsize',15);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    xlim([0 T]);
    ylabel('$c(t)$','Interpreter','latex','fontsize',15);
    title('Exposure rate','interpreter','latex','fontsize',15);
    %save the figure
    fileName3= 'Exposure rate.pdf';
    exportgraphics(gcf, fullfile(fpath, fileName3), ...
    'ContentType', 'image', 'BackgroundColor','white');


    %death due to insecticide exposure
    figure
    plot(x,DeathInsecticide,'LineWidth',LW,'Color',color(3,:),'LineStyle',':');
    xlabel('Insecticide resistance level ($x$)','Interpreter','latex','fontsize',15);
    xlim([0 1]);
    title('Insecticide exposure death','interpreter','latex','fontsize',15);
    %save the figure
    fileName4= 'Insecticide exposure death.pdf';
    exportgraphics(gcf, fullfile(fpath, fileName4), ...
    'ContentType', 'image', 'BackgroundColor','white');
    

    %survival probabilities  
    figure
    % set(gcf,'position',[100,100,600,600])  %3eCoord=largeur, 4e=hauteur

    subplot(1,2,1)
    GrNum=1;
    colormap hot
    surf(age, x, SurvProb0', 'FaceColor', 'interp', 'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    axis tight
    view(0,90);
    colorbar
    camlight left
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15),
    xlabel('Age ($a$)','Interpreter','latex','fontsize',15),
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')} ' ...
          '          Unexposed mosquitoes'], 'interpreter','latex');

    subplot(1,2,2)
    GrNum=GrNum+1;
    colormap hot
    surf(age, x, SurvProb1', 'FaceColor', 'interp', 'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud')
    axis tight
    view(0,90)
    colorbar
    camlight left
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15),
    xlabel('Age ($a$)','Interpreter','latex','fontsize',15),
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')} ' ...
          '          Exposed mosquitoes'], 'interpreter','latex');

    %save the figure
    fileName5= sprintf('C%.g_Ps%g_SurvivalProba.pdf', cm, ProbaSurviFullySensi);
    exportgraphics(gcf, fullfile(fpath, fileName5), ...
    'ContentType', 'image', 'BackgroundColor','white');
end

%% SECTION 5: Typical dynamics with a constant exposure rate C=fct(C_est_cte,cm,cval,ProbaSurviFullySensi)  
if (1)
    %setting parameters values
    C_est_cte=1;
    cm=0.6;
    cVal= zeros(1,Ny);
    ProbaSurviFullySensi= 10^-1;        
    dd0= -log(ProbaSurviFullySensi);    
    dd1= -log(1-ProbaSurviFullySensi);  
    
    %run mosquito population model
    [c,age,x,DeathInsecticide,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,...
      E0T,E1T,id0,RelativeGain,Temergence]= qIR_MainFunction(C_est_cte,cm, ...
      VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0,gamma1,mu0, ...
      mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs,cVal,k, ...
      EmergenceThreshold);

    %main figure
    qIR_MainFigure(LW,LETTERS,color,id0,Ny,T,days_per_year,t_begin_c,time, ...
      x,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,E0T,E1T, ...
      Temergence,RelativeGain);

    %save the figure
    fileName6= sprintf('FigC%.g_ProbaSurvie%g.pdf', cm, ProbaSurviFullySensi);
    exportgraphics(gcf, fullfile(fpath, fileName6), ...
    'ContentType', 'image', 'BackgroundColor','white');
end

%% SECTION 6: Typical dynamics with a variable exposure rate C=fct(C_est_cte,cm,ProbaSurviFullySensi) 
if (0)
    %setting parameters values
    C_est_cte= 0;
    cm= 0.2;               % SeuilC1<cm<SeuilC2
    ProbaSurviFullySensi=10^-10;
    dd0= -log(ProbaSurviFullySensi);    
    dd1= -log(1-ProbaSurviFullySensi);  
  
    %objective function for cVal
    objFun = @(cVal) qIR_FuncOpticVal(cVal, C_est_cte, cm, VarJ0, VarJ1, time, T, ...
                    t_begin_c, tau, age, da, x, dx, dt, x0, Na, Nx, Nt, ...
                    gamma0, gamma1, mu0, mu1, rrm, rr0, rr1, dd0, dd1, ...
                    AvrMosquiLifeSpan, MinMosquiAgeLayEggs, k, EmergenceThreshold);

    %finding Copt
    Ops=optimset('Display','iter', ...
          'Algorithm','interior-point', ...
          'MaxFunEvals',10^10, ...
          'MaxIter',20000, ...
          'TolX',0.0001, ...
          'TolFun',0.0001);
    %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

    Cinit=0.15*ones(1,10);

    tic
    Copt= fmincon(objFun, Cinit, [], [], [], [], zeros(10,1), cm*ones(10,1), [], Ops);
    toc

    % Multiple beeps to signal completion
    for i = 1:5
        beep
        pause(0.5)  % Half-second pause between beeps
    end

    cVal=Copt;
    % cVal= [0.2000  0.1894  0.1337  0.1194  0.1149  0.1169  0.1258  0.1645  0.1999  0.1000]; %cm=0.2
    % cVal= [0.1996  0.1863  0.1328  0.1166  0.1106  0.1102  0.1146  0.1255  0.1650  0.1922]; %cm=0.2 avec 'Algorithm','interior-point', 'MaxFunEvals',10^10, 'MaxIter',20000,'TolX',0.0001,'TolFun',0.0001) 
    % cVal= [0.2836  0.1487  0.1255  0.1124  0.1071  0.1069  0.1105  0.1192  0.1429  0.2855]; %cm=0.35 avec 'Algorithm','interior-point', 'MaxFunEvals',10^10, 'MaxIter',20000,'TolX',0.0001,'TolFun',0.0001)
    % cVal= [0.2843  0.1494  0.1272  0.1152  0.1114  0.1130  0.1204  0.1429  0.3007  0.2795]; %cm=0.5
    % cVal= [0.2839  0.1486  0.1255  0.1122  0.1071  0.1069  0.1106  0.1192  0.1432  0.3001]; %cm=0.5 avec 'Algorithm','interior-point', 'MaxFunEvals',10^10, 'MaxIter',20000,'TolX',0.0001,'TolFun',0.0001)

    %run mosquito population model
    [c,age,x,DeathInsecticide,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,...
      E0T,E1T,id0,RelativeGain,Temergence]= qIR_MainFunction(C_est_cte,cm, ...
      VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0, ...
      gamma1,mu0,mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan, ...
      MinMosquiAgeLayEggs,cVal,k,EmergenceThreshold);    

    % graphic representations
    LW=1.5;
    figure
    %set(gcf,'position',[100,100,800,600]);  %3eCoord=largeur, 4e=hauteur
    axes ('fontsize',15);

    %exposure rates
    ax=subplot(1,2,1);
    GrNum=1;
    plot(time,c,'LineWidth',LW,'color','k','linestyle','-')
    xlim([0 T])
    xlabel('Time $t$ (year)','Interpreter','latex','fontsize',15)
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylabel('$c(t)$','Interpreter','latex','fontsize',18)
    title(['\fontsize{13}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')} Exposure rate'], ...
        'interpreter','latex');

    %Rgain and Temg
    GrNum=GrNum+1;
    subplot(1,2,2);
    y=AA0T(id0)+AA1T(id0);
    z=AA0T+AA1T;
    plot([0 T],[y y],'LineWidth',LW,'color',color(1,:),'linestyle','-');
    xlim([0 T]);
    xlabel('Time $t$ (year)','Interpreter','latex','fontsize',15);
    xticks(days_per_year:days_per_year:T);
    xticklabels(0:1:Ny);
    ylim([min(z) max(z)]);
    ylabel('AFMs ($A_0+A_1$)','Interpreter','latex','fontsize',15)
    hold on
    plot(time,z,'LineWidth',LW,'color',color(2,:),'linestyle','--')
    TemergPlot=Temergence+t_begin_c;
    plot([TemergPlot TemergPlot],[min(z) max(z)],'LineWidth',.5,'color','k','linestyle',':')
    hold off
    title(['\fontsize{13}{0}\selectfont', '\quad', '\textbf{(' LETTERS(GrNum) ')}' ...
           '\fontsize{13}{0}\selectfont$\:{\rm r}_{\rm gain}=$' num2str(RelativeGain,2) ...
           '\ ; \fontsize{13}{0}\selectfont$\:{\rm T}_{\rm emg}=$' num2str(Temergence/300,3)], ...
           'interpreter','latex');
    legend('$c=0$','$c>0$','Interpreter','latex','location','best',...
        'Orientation','vertical','fontsize',13)
    legend boxoff 

    %save the figure
    fileName7= sprintf('FigC%.2f_NonConstant.pdf', cm);
    exportgraphics(gcf, fullfile(fpath, fileName7), ...
    'ContentType', 'image', 'BackgroundColor','white');

end
