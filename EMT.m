function [W,EE,Pmito,Hplasmy] = EMT(mu,muu,rsex,pl,plr,P,fr)
% Variables:    mu & muu - mtDNA mutation rates a->A and A->a, between 0
%               and 1
%               rsex - sex rate, between 0 and 1
%               pl - wild type paternal leakage, in multiples of 1/M,
%               between 0 and 1
%               plr - mutant paternal leakage, in multiples of 1/M, between
%               0 and 1
%               P - period of environmental fluctuations, integer
%               fr - initial frequency of the mutant, between 0 and 1

% Returns:  W - paternal leakage mutant frequeny history over time
%           EE - environmental change history over time
%           Pmito - mitochdonrial mutant distributions
%           Hplasmy - heteroplasmy distributions
%
% By default uses maternal regulation of paternal leakage. For paternal,
% see line ~136
% For mean fitness, see line ~164 and modify returns.

    M=20; % M is the number of discrete mitochondria per cell
    n=1;  % number of possible nuclear values: aa aA AA
    x=2;  % x is the exponent in the fitness function. x>1 makes it concave which is what we want.
    s=1;  % selection strength
 
    % Initialize population P
    % Maternal population, females
    P1=zeros(M+1,n);
    P1(1,1)=1;
    % Paternal population, males
    P2=zeros(M+1,n);
    P2(1,1)=(1-fr);
    % Paternal population, mutant males
    P2R=zeros(M+1,n);
    
    % Life cycle operators
    % Mitochondrial mutation
    % Both ways, rates mu and muu
    U=zeros(M+1,M+1);
    for i=0:M
        for j=0:M
            k=0:M;
            U(i+1,j+1)=sum(binopdf(k,M-j,mu).*binopdf(k+j-i,j,muu));
        end
    end
    
    % w is a matrix version of the fitness function
    % for the two environments
    w=zeros(M+1,2);
    for i=0:M
        % just count the number of mismatches
        % w = 1 - (total mismatches/M)^x;
        w(i+1,1)=1 - s*( i / M     )^x; % environment a (1)
        w(i+1,2)=1 - s*( (M-i) / M )^x; % environment A (2)
    end
    
    % Resampling at Meiosis 1
    L1=zeros(2*M+1,2*M+1);
    for i=0:2*M
        for j=0:2*M
            L1(i+1,j+1)=hygepdf(i,4*M,2*j,2*M);
        end
    end
    
    % Resampling at Meiosis 2
    L2=zeros(M+1,2*M+1);
    for i=0:M
        for j=0:2*M
            L2(i+1,j+1)=hygepdf(i,2*M,j,M);
        end
    end
    
    % Mitochondrial resampling after cell fusion
    % With paternal leakage pl, the zygote has M+pl*M total mitochondria
    % PSI resamples that up to 2M, a proper number for the zygote
    PSI=zeros(0,1);
    for i=0:2*M
        for j=0:round(pl*M)+M
            PSI(i+1,j+1)=binopdf(i,2*M,j/(round(pl*M)+M));
        end
    end
    
    % Mitochondrial resampling at cell fusion for the paternal cell
    % With paternal leakage it should only contribute pl*M mitochondria
    % PSIP resamples the number of mito mutants without replacement
    PSIP=zeros(0,1);
    for i=0:pl*M
        for j=0:M
            PSIP(i+1,j+1)=hygepdf(i,M,j,round(pl*M));
        end
    end
    
    PSIR=zeros(0,1);
    for i=0:2*M
        for j=0:round(plr*M)+M
            PSIR(i+1,j+1)=binopdf(i,2*M,j/(round(plr*M)+M));
        end
    end
    PSIPR=zeros(0,1);
    for i=0:plr*M
        for j=0:M
            PSIPR(i+1,j+1)=hygepdf(i,M,j,round(plr*M));
        end
    end
    
    % Mitochondrial resampling for asexual doubling and clonal division
    A=zeros(0,1);
    for i=0:M
        for j=0:M
            A(i+1,j+1)=hygepdf(i,2*M,2*j,M);
        end
    end

    E = 0;
    % Now iterate the life cycle for many many generations k
    for k=1:(50100)
        if k==1
            P2R=fr*P2;
            P2=P2-fr*P2;
        end
        if rand()<1/P
        %if k==5000
            E=abs(E-1);
        end
        % Mutate mitochondria from a to A, from A to a
        P1 = U*P1;
        P2 = U*P2;
        P2R = U*P2R;
       
        % Selection
        ww = w(:,E+1);
        P1 = (ww.*P1) / sum(ww.*P1,'All');
        P2P2R = ([ww,ww].*[P2,P2R]) / sum([ww,ww].*[P2,P2R],'All');
        P2 = P2P2R(:,1);
        P2R = P2P2R(:,2);
        
        % Sexual distributions S
        % maternal regulation of mt inheritance
        Z=PSI*(conv(PSIP*P1,P2)*2);
        ZR=PSIR*(conv(PSIPR*P1,P2R)*2);
        % paternal regulation of mt inheritance
        %Z=PSI*(conv(P1,PSIP*P2)*2);
        %ZR=PSIR*(conv(P1,PSIPR*P2R)*2);

        
        Zsum = sum(Z+ZR,'All');
        Z = Z/Zsum;
        ZR = ZR/Zsum;       

        S1=L2*L1*( Z+ZR );   
        S2=L2*L1*( Z );
        S2R=L2*L1*( ZR );

        % asexual and sexual combined
        P1=A*P1*(1-rsex)+S1*rsex;
        P2=A*P2*(1-rsex)+S2*rsex;
        P2R=A*P2R*(1-rsex)+S2R*rsex;
   
        % if we need to keep mutant frequency constant to measure fitness
        % advantage, set it to fr here:
        %P2 = P2/sum(P2,'All')*(1-fr);
        %P2R = P2R/sum(P2R,'All')*(fr);

        % record mean fitness or frequencies
        W1(k) = sum(ww.*P1,'All');
        W2(k) = sum([ww,ww].*[P2,P2R],'All');
        
        %W(k) = sum(ww.*P2R,'All')/sum(P2R,'All')-sum(ww.*P2,'All')/sum(P2,'All') ;
        W(k) = sum(P2R,'All');
        EE(k) = E;

    end

    % mitochondrial distributions
    Pmito = [sum(P1,2)/sum(P1,'all'),sum(P2,2)/sum(P2,'all'),sum(P2R,2)/sum(P2R,'all')];
    H1 = [P1(1:10)+flip(P1(12:end));P1(11)];
    H2 = [P2(1:10)+flip(P2(12:end));P2(11)];
    H3 = [P2R(1:10)+flip(P2R(12:end));P2R(11)];
    Hplasmy = [H1/sum(H1,'all'),H2/sum(H2,'all'),H3/sum(H3,'all')];
    
    %S=var(W(end/5:end));
    %W=mean(W(end/5:end));
    
end
    