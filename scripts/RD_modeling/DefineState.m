function [y,t,s] = DefineState(E)
%This is a function that uses the eigenvalues of the Jacobian evaluated
%without diffusion to define which system type (0,1,2,3) is present.

%Input(s):
% - E: eigenvalues of Jacobian calculated without diffusion

%Output(s):
% - y: Binary value that defines whether a system should be further
% analysed for Turing Instabilities
% - t: threshold value
% - s: system type (0 = stable, 1 = unstable, 2 = oscillator, 3 = damped
% oscillator

    if all(real(E) < 0) && isreal(E) == 1 %This is stable
        s = 0;
        y = 1;
        t = 0;
    elseif all(real(E) < 0) && isreal(E) == 0 %This is a damped oscillator
        s = 3;
        y = 1;
        t = 0;
    elseif all(real(E) == 0) && isreal(E) == 0 %This is a limit cycle
        s = 2;
        y = 0;
        t = max(real(E(:)));
    else
        s = 1; %This is unstable
        y = 0;
        t = 0;
    end
end