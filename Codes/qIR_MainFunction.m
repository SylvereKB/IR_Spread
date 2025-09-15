function[c,age,x,DeathInsecticide,Theta0,Theta1,AA0,AA1,E0,E1,AA0T,AA1T,...
    E0T,E1T,id0,RelativeGain,Temergence]= qIR_MainFunction(C_est_cte,cm, ...
    VarJ0,VarJ1,time,T,t_begin_c,tau,age,da,x,dx,dt,x0,Na,Nx,Nt,gamma0,gamma1,mu0, ...
    mu1,rrm,rr0,rr1,dd0,dd1,AvrMosquiLifeSpan,MinMosquiAgeLayEggs,cVal,k, ...
    EmergenceThreshold)
        
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
