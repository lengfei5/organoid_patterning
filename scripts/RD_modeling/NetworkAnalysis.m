%Run SystemDefinition_example first

%%
%__________________________________________________________________________
%ODE optimisation parameters
options = optimoptions('fsolve','FunctionTolerance', 1.0e-12,'Display','off'); %ODE solver setting see, Matlab documentation for further details
t_final = 1000; %Time threshold for ODE simulation
%__________________________________________________________________________
evaluate_diffusion = 1; %Set to 1 if Jacobian should be evaluated with diffusion. 
% If set to 0 only the Jacobian excluding diffusion will be evaluated and no Turing instabilities will be analysed.

%%
%__________________________________________________________________________
%Initialise some matrices
y0 = rand.*ones(n, 1); %Initial conditions for first simulation     
total = length(k_grid(1,:)); %Calculates number of iterations needed to sample all parameter combinations
z = k_length+2; 
System_save = zeros(total,z); %System save is a matrix to save all parameters as well as information to the amount of Steady States 
%found for each parameter set and the classfication of the system according to 0(Stable), 2(unstable), 3(oscillating), 4(damped oscillating). 
%The behavior is estimated from the simulation trajectories.
SteadyState_save = zeros(2*total,n+3); %This is a matrix that saves all steady states indexed to each parameter set.
Turing_save = zeros(2*total,n+3); %A matrix to save which parameter combination and corresponding steady state form Turing Instabilities. 
% In addition, the diffusion values, and what type of Turing Instability is present, are defined.
tspan = [0 t_final]; %Vector defining time for solving ODE equation using the Matlab ODE suite

count_i = 1;
len = length(c_ini(:,1));
Saver=zeros(len^n,n);
counter = 2+n+k_length; 
final = counter;

%%
%__________________________________________________________________________
%Workflow Part 1: Find steady states numerically
tic;

for i = 1:length(k_grid(1,:)) % loop over the sampled parameters 
    k = k_grid(:,i)'; %Select k's to be sampled from k_grid  
    try [t1, xout1] = ode15s(@(t, x)ode(t,x,f_ode,k),tspan,y0); %Run first ODE simulation
    catch
        [t1, xout1] = ode23s(@(t, x)ode(t,x,f_ode,k),tspan,y0);
    end
    
    % figure; plot(t1,xout1(:, 1), t1, xout1(:, 2))
    ss2 = xout1(length(xout1(:,1)),:); %Suggest steady state solution from simulation
    m = BurnInCorrection(t1,xout1); %Calculate threshold for local maxima estimation
    
    state = 0;
    
    %Check if system is showing any oscillations:
    max_min_count = 0;
    for p1 = 1:n
        local_max = (findpeaks(xout1(:,p1),'MinPeakProminence',m(p1)));
        if length(local_max) > 1
            max_min_count = max_min_count + 1;
            local_max_end(p1) = local_max(length(local_max)-1);
            local_max_save = local_max;
        end 
    end
    
    %Distinguish damped oscillator from oscillator using the ODE trajectories
    if max_min_count > 0
        state = OscVsDamped(local_max_save(length(local_max_save)-1),local_max_save(length(local_max_save)));
        ss = median(xout1(:,:));
    end
    %Evaluate if steady state guess is correct for non-oscillating systems
    
    %Check for negative and imaginary solution and exclude
    ss2(ss2<=0)=0;
    if all(ss2) == 0 || isreal(ss2) < 1
      state = 1; 
    end
    
    if state == 0 %For systems with no oscillatory behaviour and that have been found to be stable.
        %Checking for multiple steady states simulating several ODE
        %trajectories for all initial conditions defined in c_ini
        Saver = ODEmulti(c_ini,n,k,f_ode,tspan); %ODEmulti performs the ODE simulations and exports all endpoints of the trajectories
        %Perform cluster analysis to check for redundant entries
        C = ClusterAnalysis(Saver,n,k,f_ode,tspan);
        if isempty(C) == 1
            state = 1;
        end
        
        %Save all non-redundant steady states
        for p1 = 1:length(C(:,1))
            System_save(i,length(k)+2)=length(C(:,1)); %save number of steady states (e.g. 2 = 2 steady states)
            SteadyState_save(count_i,1) = i;
            SteadyState_save(count_i,2:n+1)=C(p1,:); %save steady state (steady state values for each node)
            count_i = count_i + 1;
        end
        
    end

    if state == 2 %System is found to be an oscillator and will not be further analysed for any Turing Instabilities
        SteadyState_save(count_i,1) = i;
        SteadyState_save(count_i,2:n+1)=ss;
        SteadyState_save(count_i,n+2)=state;
        count_i = count_i + 1;
    end
    
    if state == 3  %System is defined a damped oscillating system
        %Checking for multiple steady states by simulating all trajectories
        %defined in c_ini
        Saver = ODEmultiOsc(c_ini,n,k,f_ode,tspan);
        %Perform cluster analysis: using f_solve, solutions are optimised
        %and subsequently k_means clustering is performed to distinguish
        %steady states.
        C = ClusterAnalysisOsc(Saver,n,k,f_ode,options);
        System_save(i,length(k)+2)=length(C(:,1)); %save number of steady states
        for p1 = 1:length(C(:,1))        
            SteadyState_save(count_i,1) = i;
            SteadyState_save(count_i,2:n+1)=C(p1,:); %save steady state values for each ndoe
            count_i = count_i + 1;
        end   
    end
        
    if state ==1 %if systems does not converge sufficienlty within the time given in tspan, define system as unstable.
        System_save(i,length(k)+2)=0;
    end
    
    System_save(i,1:length(k))=k;
    System_save(i,length(k)+1)=state;
end

%Get rid of all empty entries
for i = 1:length(SteadyState_save(:,1))
    if SteadyState_save(i,1) ==0 
        j = i-1;
        break
    end
end

SteadyState_save = SteadyState_save(1:j,:);

toc

%% Workflow Part II:
%__________________________________________________________________________
% Given all steady states that are either stable or damped 
% oscillators the Jacobian eigenvalues are calculated without diffusion. 
%For all systems where maximum(real(eigenvalue)) < 0 the eigenvalues for the
%Jacobian including diffusion are then analysed. If any
%maximum(real(eigenvalue)) > 0 the system exhibits a Turing I instability.
%Depending on the behaviour of the instability (finite or infinite qmax) a
%Turing I or II instability is classified respectively. 

tic %Timing process 2


% here is the symbolic part making trouble
X = sym('x', [1 n]); %Create symbolic nodes x1...xn
K = sym('k',[1 k_length]); %Create symbolic parameters k1...kn
f_sym = f_ode(X,K); %Create symbolic function of f_ode
Z = jacobian(f_sym,X); %Calculate Jacobian matrix
func_Z = matlabFunction(Z); %Transfer symbolic Jacobian to matlab function
Zinputs=symvar(Z); % read out input parameters of jacobian Z
Xoverlap = ismember(X,Zinputs); %find nodes that have to be substituted in Jacobian
Koverlap = ismember(K,Zinputs); %find parameters that have to be substituted in Jacobian

count = 0;

for i = 1:length(SteadyState_save(:,1))
    if System_save(SteadyState_save(i,1),length(k)+1) ~= 1 && SteadyState_save(i,n+2) ~= 2 %Exclude oscillators and unstable solutions
        kk = System_save(SteadyState_save(i,1),1:k_length);
        xx = SteadyState_save(i,2:n+1);
        input = num2cell([kk(Koverlap), xx(Xoverlap)]); 
        S = func_Z(input{:}); %Substitute steady state values for nodes and parameters in Jacobian
        E = eig(S); %Calculate eigenvalues WITHOUT diffusion
        [y,t,s] = DefineState(E); %Test if all real eigenvalues are negative and whether complex part exists in eigenavlue, define state accordingly
        %Correct Steady State count when system unstable defined by
        %max(real(E))> 0
        if s == 1 && System_save(SteadyState_save(i,1),k_length+2) ~= 0 %if egenvalues define system as unstable --> delete entry from steady state entries.
            System_save(SteadyState_save(i,1),k_length+2) = System_save(SteadyState_save(i,1),k_length+2)-1;
        end
        SteadyState_save(i,n+2) = s;
        SteadyState_save(i,n+3) = max(real(E));
        if evaluate_diffusion == 1 && y == 1 %If diffusion should be evaluated (default setting) AND the system converges to a steady state.
            for j = 1:length(binary_diffusor(:,1))
                diffusing_nodes = binary_diffusor(j,:); %Call in definition of which node diffuses from matrix/vector binary_diffusor
                [w1,D_matrix] = CreateDiffusionMatrix(n,m_min,m_max,d_range,diffusing_nodes); %This function will output a diffusion matrix if the node number of diffusing nodes is within the range of m_min and m_max. In addition a binary output w1 is parsed to the following if loop.
                %Analysis for Turing instabilities:
                if w1 == 1
                    for p1 = 1:length(D_matrix(1,:))
                        D = D_matrix(:,p1);
                        w2 = IsTuring(S,D,t); %Function IsTuring takes in the substituted form of the Jacobian matirx as well as a vector defining the diffusion parameters for each node. Next eigenvalues are calculated of the Jacobian INCLUDING diffusion. w2 defines if and which type of Turing instability (1 or 2) has been found.
                        %Save what parameter combination, steady state and
                        %diffusion values exhibit a Turing I instability.
                        if w2 == 1
                            count = count + 1;
                            Turing_save(count,1) = i; %Index corresponds to steady state
                            Turing_save(count,2:n+1) = D; %Save diffusion values
                            Turing_save(count,n+2)= j; %save "type" of diffusors = which are the diffusing nodes (indexed to binary diffusor vector/matrix)
                            Turing_save(count,n+3)= w2; %save Type of Turing Instability                   
                        elseif w2 == 2
                            count = count + 1;
                            Turing_save(count,1) = i; %Index corresponds to steady state
                            Turing_save(count,2:n+1) = D; %Save diffusion values
                            Turing_save(count,n+2)= j; %save "type" of diffusors = which are the diffusing nodes (indexed to binary diffusor vector/matrix)
                            Turing_save(count,n+3)= w2; %save Type of Turing Instability  
                        end                     
                    end
                end
            end
        end
    end
end

toc

%% Delete all zero entries and save 
for i = 1:length(Turing_save(:,1))
    if Turing_save(i,1) ==0 
        j = i-1;
        break
    end
end

Turing_save = Turing_save(1:j,:); %Get rid of unused entries except for 1 

%Save outputs to a file named ID (see Systemdefinition_example.m)
save(ID,'System_save','SteadyState_save','Turing_save')

