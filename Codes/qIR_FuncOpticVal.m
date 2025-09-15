function[RelativeGain2,Temergence]=qIR_FuncOpticVal(cVal,C_est_cte,cm, ...
      VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0, ...
      gamma1,mu0,mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan, ...
      MinMosquiAgeLayEggs,k,EmergenceThreshold)

    %% This Function computes the relative gain and the time of emergence of IR
    %
    % Inputs:
    %   cVal - vector of insecticide exposure per year
    %   C_est_cte - boolean value 1 if exposure rate is constant over deploymenbt period, 0 if not
    %   cm - insecticide exposure rate
    %   VarJ0 - mutational variance of reference susceptible mosquito
    %   VarJ1 - mutational variance of reference resistant mosquito
    %   time - vector of discretisation of time
    %   T - total time
    %   t_begin_c - starting time of insecticide exposure
    %   tau - proportion of hatched eggs laid by AFM that reach adulthood
    %   age - discretisation of AFM age
    %   da - age step
    %   x - discretisation of insect. resis. level 
    %   dx - insect. resis. step
    %   dt - time step
    %   x0 - reference sensitive  insecticide resistance level
    %   Na - number of age steps
    %   Nx - number of insect. resist. level steps
    %   Nt - number of time steps
    %   gamma_i - hatching rate of eggs laid by unexposed/exposed AFM
    %   mu_i - natural death rate of eggs laid by unexposed/exposed AFM
    %   rrm - maximum number of eggs laid by the ref sensi strain
    %   rr0 - average number of eggs laid by the ref sensi strain
    %   rr1 - average number of eggs laid by the ref resis strain
    %   dd0 - death rate of ref sensi strain due to insecticide
    %   dd1 - death rate of ref resis strain due to insecticide
    %   AvrMosquiLifeSpan - average mosquito life span in days
    %   MinMosquiAgeLayEggs - minimum age for AFM to lay eggs in days
    %   k - exposant of the logistic function growth of eggs 
    %   EmergenceThreshold - threshold in % of detection of 
    %
    % Outputs:
    %   RelativeGain2 - relative gain of using insecticide durint the deplyment period 
    %   Temergence - time at which 10% of mosquito are now resistant  

    %% Run qIR_MainFunction
    [~,~,~,~,~,~,~,~,~,~,AA0T,AA1T,~,~,id0,~,~]= qIR_MainFunction(C_est_cte,cm, ...
    VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0,gamma1,mu0, ...
    mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs,cVal,k, ...
    EmergenceThreshold);

    %% Computation of Temergence
    resistance_percentage = 100*AA1T./(AA0T+AA1T);

    [Row_Temg,~]= find(resistance_percentage > EmergenceThreshold);

    if isempty(Row_Temg)  
       id_Temg= Nt;
    else 
        GRAD= gradient(Row_Temg);
        id_Temg= Row_Temg(1);
        for np=2:length(Row_Temg)
            if GRAD(np)>1, id_Temg= Row_Temg(np); 
            end
        end
    end

    Temergence= time(id_Temg)-t_begin_c;
    
    %% Computation of Relative Gain
    TimeDdpl= t_begin_c:dt:T;

    IdTimeDdpl= (1+floor(t_begin_c/dt)):Nt;

    RelativeGain2= trapz(TimeDdpl,AA0T(IdTimeDdpl)+AA1T(IdTimeDdpl))/...
        ((AA0T(id0)+AA1T(id0))*(T-t_begin_c));

end
