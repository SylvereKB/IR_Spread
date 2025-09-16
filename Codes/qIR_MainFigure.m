function[c,age,x,DeathInsecticide,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,...
    E0T,E1T,id0,RelativeGain,Temergence]= qIR_MainFunction(C_est_cte,cm, ...
    VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0,gamma1,mu0, ...
    mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs,cVal,k, ...
    EmergenceThreshold)
        
    %% This Function gives the total number of all subpop of mosquitoes, the relative gain and the time of emergence of IR
    %
    % Inputs:
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
    %   cVal - vector of insecticide exposure per year
    %   k - exposant of the logistic function growth of eggs 
    %   EmergenceThreshold - threshold in % of detection of IR
    %
    % Outputs:
    %   c - insecticide exposure rate in days
    %   age - vector of discretisation of age
    %   x - vector of discretisation of insec. resis. level
    %   DeathInsecticide - death due to insecticide exposure
    %   Theta0 - fitness function of ref sensitive mosquito population
    %   Theta1 - fitness function of ref resitant mosquito population
    %   AA0 - number of unexposed AFM in function of time and insec. resis. level
    %   AA1 - number of exposed AFM in function of time and insec. resis. level
    %   E0 - number of unexposed eggs in function of time and insec. resis. level
    %   E1 - number of exposed eggs in function of time and insec. resis. level
    %   AA0T - number of unexposed AFM in function of time
    %   AA1T - number of exposed AFM in function of time
    %   E0T - number of unexposed eggs in function of time
    %   E1T - number of exposed eggs in function of time
    %   id0 - indice of t_begin_c
    %   RelativeGain - 'normalized' relative gain of using insecticide durint the deplyment period 
    %   Temergence - time at which 10% of mosquito are now resistant  

    %% Mutation kernels of unexposed/exposed AFM
    FunJ0= @(x)normpdf(x,0,VarJ0);     %mutation function from allele x to y 
    FunJ1= @(x)normpdf(x,0,VarJ1);     %mutation function from allele x to y 
    J0= FunJ0(dx*(1-Nx:Nx-1));         %mutation kernel from allele x to y
    J1= FunJ1(dx*(1-Nx:Nx-1));         %mutation kernel from allele x to y
    
    %% Insecticide exposure rate
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
    
    %% Computation of r0, r1 and Theta0 and Theta1 
    [Theta0,Theta1,r0,r1,~,~,~,~,DeathInsecticide,d0,d1]=...
      qIR_FuncTheta(tau,age,x,Na,Nx,gamma0,gamma1,mu0,mu1,rrm,rr0,rr1,dd0, ...
      dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs);    
    
    % %% Computation of spectral radii of operators L0, L1 and L
    % R_0 = rayon_spectral(FunJ0,Theta0,x,dx); %=1.7243e+03
    % R_1 = rayon_spectral(FunJ1,Theta1,x,dx); %=1.5539e+03
    % SeuilC1=1-R_1/R_0;%= 0.0988
    % SeuilC2=1-1/R_0;%= 0.9994
    
    %% Computation of E0, E1, A0 and A1.
    %Initialization 
    E0= zeros(Nt,Nx);                 A0= zeros(Na,Nx,Nt);
    E1= zeros(Nt,Nx);                 A1= zeros(Na,Nx,Nt);

    E0(1,:)= 1*normpdf(x,x0,0.01);
    
    for t=1:(Nt-1)
        %computation of H(E(t))
        A0A1Int=zeros(1,Nx);
        for idx=1:Nx
            A0A1Int(idx)=trapz(age,(A0(:,idx,t)+A1(:,idx,t))/tau(idx));
        end
        Et=trapz(x,E0(t,:)+E1(t,:)+A0A1Int);
        HEt=(1+Et)^(-k);
    
        %computation of \int_Omega m0(x-y)[\int_0^\infty r0(a,y)*A0(a,y,t)da]dy
        %and \int_Omega m1(x-y)[\int_0^\infty r1(a,y)*A1(a,y,t)da]dy
        r0A0= trapz(age,r0.*A0(:,:,t));
        r1A1= trapz(age,r1.*A1(:,:,t));
        
        ConvoKernel0=dx*conv(r0A0,J0,'same');  
        ConvoKernel1=dx*conv(r1A1,J1,'same');  
    
        %computation of E0 and E1
        E0(t+1,:)= (E0(t,:)/dt + HEt*ConvoKernel0)./(1/dt + mu0 + gamma0);
        E1(t+1,:)= (E1(t,:)/dt + HEt*ConvoKernel1)./(1/dt + mu1 + gamma1);
    
        %computation of A0 and A1
        A0(1,:,t+1)= (1-c(t))*tau.*gamma0.*E0(t+1,:);
        A1(1,:,t+1)= c(t)*tau.*gamma0.*E0(t+1,:) + tau.*gamma1.*E1(t+1,:);

        for a=2:Na
            A0(a,:,t+1)=(A0(a,:,t)/dt + A0(a-1,:,t+1)/da)./...
                (1/dt + 1/da + d0(a,:));
            A1(a,:,t+1)=(A1(a,:,t)/dt + A1(a-1,:,t+1)/da)./...
                (1/dt + 1/da + d1(a,:));
        end     
    end
    
    
    %% Computation of total pop of eggs and AFM at time t
    %computation of total eggs population evolution    
    E0T= squeeze(trapz(x,E0,2));
    E1T= squeeze(trapz(x,E1,2));

    %computation of total AFM population evolution
    %size A0=(a,x,t). trapz(age,A0,1) return a matrix of size (1,x,t).
    %squeeze(trapz(age,A0,1)) return a matrix of size (x,t) ie supression of first dim
    %squeeze(trapz(age,A0,1))' return the transpose matrix of size (t,x)
    AA0= (squeeze(trapz(age,A0,1)))';
    AA1= (squeeze(trapz(age,A1,1)))';
       
    AA0T= squeeze(trapz(x,AA0,2));
    AA1T= squeeze(trapz(x,AA1,2));

    %% Computation of Temergence
    resistance_percentage = 100*AA1T./(AA0T+AA1T);
    [Row_Temg,~]= find(resistance_percentage>EmergenceThreshold);
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

      % % Add this after the Temergence calculation
      % figure;
      % plot(time, resistance_percentage);
      % hold on;
      % yline(EmergenceThreshold, 'r--', 'Threshold');
      % xline(time(id_Temg), 'k--', 'Detected Emergence');
      % xlabel('Time'); ylabel('Resistance %');
      % title('Resistance Evolution');

    %% Computation of Relative Gain
    TimeDdpl= t_begin_c:dt:T;
    IdTimeDdpl= (1+floor(t_begin_c/dt)):Nt;
    RelativeGain= 1-trapz(TimeDdpl,AA0T(IdTimeDdpl)+AA1T(IdTimeDdpl))/...
        ((AA0T(id0)+AA1T(id0))*(T-t_begin_c));

end
