function w = IsTuring(S,D,threshold)
% This is a function that tests if a system exhibits a Turing Instability
% by calculating the behaviour in respect to wavenumber q of the eigenvalues 
% for the Jacobian matrix evaluated with diffusion.

% Output(s): 
% - w = 0 for no Turing instability
% - w = 1 for a Turing I instability
% - w = 2 for a Turing II instability

% Input(s): 
% - S (Jacobian evaluated at Steady State
% - D (Vector defining diffusion rates of nodes)
% - threshold (should be set to 0)

    q = logspace(-4,8,50); %Define wavenumbers to sample
    for j = 1:50
       S2=S-diag(D*(q(j))); %Calculate Jacobian matrix including diffusion
       E = eig(S2); %calculate eigenvalues
       Eig_save(j) = max(real(E)); %find maximum eigenvalue
    end
    if max(Eig_save) > threshold %Check if any positive real eigenwert exists 
        if Eig_save(length(Eig_save)) < max(Eig_save)-0.01*max(Eig_save)
            %Check if Type I
            w = 1;
        else
            %This is Type II
            w = 2;
        end
    else 
        %No Turing instability
        w = 0;
    end
    return
end