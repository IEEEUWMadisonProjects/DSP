function [ phaseShift, amplitude ] = sineFit( s, t, omega )
%PHASESHIFT This function is used to calculate the phase differences
%between four different sinusoidal inputs.
%   The chosen model of the sinusoid to fit is y = D + B*sin(omega*t+phi)
%   Using the in-phase and quadrature basis, we know that y can be
%   represented like y = D + b*cos(omega*t) + a*sin(omega*t). Omega is
%   assumed to be known. phi = atan(a/b) B = 2norm(a,b). We are looking to
%   find D, b, and a.
% Author: Alex Gabourie
%initialize A matrix of Ax = s 
A = zeros(length(s), 3);
A(:,1) = 1;
A(:,2) = cos(omega*t);
A(:,3) = sin(omega*t);

% We now determine the phase of the signal
A_n = A'*A;
x1 = A_n\(A'*s');
%returns phase between -pi and pi
phaseShift = atan2(x1(2),x1(3));
amplitude = sqrt(x1(2)^2+x1(3)^2);

end
