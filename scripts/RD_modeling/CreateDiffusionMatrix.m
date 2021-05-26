function [w, D_matrix] = CreateDiffusionMatrix(n,m_min,m_max,d_range,diffusing_nodes)
%This is a function that creates a matrix of diffusion parameters to be
%sampled. If the amount of diffusing_nodes is set within the defined range 
%of m_min and m_max this combination of diffusing nodes is allowed. 
%Using d_range all entries to be sampled are calculated. 
%Every first node that diffuses (e.g. for diffusing_nodes = [0 1 1] the 
%first node to diffuse is node 2) is set to a diffusion value of 1, 
%all other nodes are then sampled at ranges defined in d_range. 

%Input(s):
% - n = node number
% - m_min = minimum number of diffusing nodes
% - m_max = maximum number of diffusing nodes
% - d_range = range of diffusion values that are to be sampled (for
% diffusing nodes 2...m
% - diffusing_nodes: binary vector to define which nodes
% diffuse (= 1's) or do not diffuse (= 0's)

%Output(s)
% - w = binary value that is set 0 if no diffusing
% - D_matrix = minimum number of diffusing nodes

a_0 = [0]; %For nodes that do not diffuse, diffusion value is set to 0
a_1 = [1]; %For first node to diffes the diffusion value is set to 1
a_m = d_range; %For all non first-nodes to diffuse d is sampled in range d_range
for jj = m_min:m_max
    if sum(diffusing_nodes) == jj
        count = 0;
        for jjj = 1:length(diffusing_nodes)
           if diffusing_nodes(jjj) == 0 %non diffusing nodes
               d{jjj} = a_0;
           elseif diffusing_nodes(jjj) == 1 && count == 0 %first diffusing node
               d{jjj} = a_1;
               count = 1;
           else %2...nth diffusing node
               d{jjj} = a_m;
           end
        end
    end
end
if sum(diffusing_nodes) >= m_min && sum(diffusing_nodes) <= m_max
    D_matrix = combvec(d{:}); %Create combinatorial matrix from all diffusion values
    w = 1; %Set w = 1 so that diffusion is evaluted in main script
else
    D_matrix = [0]; %Diffusing nodes not within range m_min : m_max so D_matrix is empty;
    w = 0; %Return w = 0, no diffusion has to be evaluted in main script
end

end