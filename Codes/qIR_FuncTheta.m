function[Theta0,Theta1,r0,r1,r0X,r1X,SurvProb0,SurvProb1,DeathInsecticide,d0,d1]=...
   qIR_FuncTheta(tau,age,x,Na,Nx,gamma0,gamma1,mu0,mu1,rrm,rr0,rr1,dd0,dd1, ...
   AvrMosquiLifeSpan,MinMosquiAgeLayEggs)

    %% This Function computes fitness functions and Survival probabilities
    %
    % Inputs:
    %   tau - proportion of hatched eggs laid by AFM that reach adulthood
    %   age - discretisation of AFM age
    %   x - discretisation of insect. resis. level 
    %   Na - number of age steps
    %   Nx - number of insect. resist. level steps
    %   gamma_i - hatching rate of eggs laid by unexposed/exposed AFM
    %   mu_i - natural death rate of eggs laid by unexposed/exposed AFM
    %   rrm - maximum number of eggs laid by the ref sensi strain
    %   rr0 - average number of eggs laid by the ref sensi strain
    %   rr1 - average number of eggs laid by the ref resis strain
    %   dd0 - death rate of ref sensi strain due to insecticide
    %   dd1 - death rate of ref resis strain due to insecticide
    %
    % Outputs:
    %      Theta_i - Fitness function of unexposed/exposed Adult Mosquitoes
    %            r - laying rate of unexposed/exposed AFM
    %           rX - egg-laying rate as fct of insecticide resistance level
    %  Survproba_i - mortality rate of unexposed/exposed AFM
    %  DeathInsecticide - death rate as fct of IR level
    %          d_i - mortality rate of unexposed/exposed AFM
    
    %% Mortality rate of unexposed/exposed AFM
    Bard0= 10*(age>AvrMosquiLifeSpan);  
    d0= zeros(Na,Nx);
    d1= zeros(Na,Nx);
    DeathInsecticide= zeros(1,Nx);
    for idx=1:Nx
        DeathInsecticide(idx)= dd0*(dd1/dd0)^(x(idx));
        d0(:,idx)= Bard0(:);
        d1(:,idx)= Bard0(:) + DeathInsecticide(idx);
    end 

    %% Computation of survival probability of unexposed/exposed AFM
    %calculating cumulative mortality integrals
    CumD0= cumtrapz(age,d0);      % Na x Nx
    CumD1= cumtrapz(age,d1);

    %survival probability = exp(-∫ mortalité)
    SurvProb0= exp(-CumD0);        % Na x Nx
    SurvProb1= exp(-CumD1);

    %% Computation of laying rate of unexposed/exposed AFM
    Barr= age>MinMosquiAgeLayEggs;
    r0= zeros(Na,Nx);             r1= zeros(Na,Nx);
    for s=1:Nx
        r0(:,s)= rrm/(1+(rrm/rr0-1)*...
            (rr0*(rrm-rr1)/(rr1*(rrm-rr0)))^x(s))*Barr(:);
    end  
    r1= r0;

    %% Egg-laying rate as function of insecticide resistance level
    r0X=squeeze(trapz(age,r0,1));
    r1X=squeeze(trapz(age,r1,1));
    
    %% Computation of Theta0 and Theta1
    Theta0= tau.*gamma0./(mu0+gamma0).*squeeze(trapz(age,r0.*SurvProb0,1)); 
    Theta1= tau.*gamma1./(mu1+gamma1).*squeeze(trapz(age,r1.*SurvProb1,1)); 

end
