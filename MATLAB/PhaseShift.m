function [ phaseShift ] = PhaseShift( s1, t, omega )
%PHASESHIFT This function is used to calculate the phase differences
%between four different sinusoidal inputs.
%   The chosen model of the sinusoid to fit is y = D + B*sin(omega*t+phi)
%   Using the in-phase and quadrature basis, we know that y can be
%   represented like y = D + b*cos(omega*t) + a*sin(omega*t). Omega is
%   assumed to be known. phi = atan(a/b) B = 2norm(a,b). We are looking to
%   find D, b, and a.

%initialize A matrix of Ax_i = s_i, where i is 1,2 
A = zeros(length(s1), 3);
A(:,1) = 1;
A(:,2) = cos(omega*t);
A(:,3) = sin(omega*t);

% We now determine the phases of the two signals
A_n = A'*A;
x1 = A_n\(A'*s1');
%x2 = A_n\(A'*s2');

% The input s1 will be the reference signal and the output phase shift will
% be the phase of s2-s1.
phaseShift = atan(x1(2)/x1(3));% - atan(x2(2)/x2(3));

end
