%%%First define the system of interest:
%Define name for saving data:       
ID=strcat('/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_test/2N_testExample_',datestr(now,30)); %name of file that will be saved.

n=2; %n specifies the node number
k_length = 10; %Number of parameters to be sampled

%% define reactions of ODE 
%Specify your ODE function:
f_ode=@(x,k)[k(1).*k(2).^(-2).*x(1).^2.*(1+k(2).^(-2).*x(1).^2+k(3).^(-2).*x(2).^2).^(-1)+ k(6)+(-1).*k(8).*x(1); 
             k(4).*k(5).^(-2).*x(1).^2.*(1+k(5).^(-2).*x(1).^2).^(-1)+ k(7) + (-1).*k(9).*x(2)];         

 %__________________________________________________________________________
%Define parameter space for sampling as k_grid
%Example:
ks = logspace(-1,2,3); %define the range and interval parameters should be sampled at

for i = 1:k_length
    k_input{i} = ks;
end
k_grid = combvec(k_input{:}); %Create a combinatorial matrix from the above defined ranges.

%Please note: k_grid should not exceed 500000 different combinations.

%%
%__________________________________________________________________________
%Define diffusion parameters: the algorithm requires four inputs for this
m_min = 2; %m_min specifies the minimum number of nodes that can diffuse
m_max = n; %m_max specifies the maximum number of nodes that can diffuse

%Here we take all possible permutations of nodes HOWEVER only systems that fulfil the m_min and m_max requirement will be sampled.
% binary_diffusor = permn([0 1],n); %Define which nodes diffuse (1 specifying diffusing, 0 non-diffusing). 
%If only specific nodes should diffuse set binary_diffurso to e.g. [1 1 0] for a three node system (A B C) where A and B diffuse. 

binary_diffusor = [1, 1]; % alternative way to specify the diffusible nodes with a vector
d_range = logspace(-3,2,5); %This is the diffusion range that will be analysed.

%Note: always adjust m_min and m_max according to the node number that
%should diffuse even if binary_diffuser consists of only one element. 

%%
%__________________________________________________________________________
%Define initial conditions for ODE
x_max = 20001;
int = 10000;
c_ini = permn(1:int:x_max,n); %Grid from which algorithm will sample. Here 3 initial conditions per node are sampled.





