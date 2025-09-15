function[RelativeGain2,Temergence]=qIR_FuncOpticVal(cVal,C_est_cte,cm, ...
      VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0, ...
      gamma1,mu0,mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan, ...
      MinMosquiAgeLayEggs,k,EmergenceThreshold)

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
