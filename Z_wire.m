% Impedancemodell with Besselfun. (function)
function [z_math] = Z_wire (omega, R, sigma, u0, ur)
  A = sqrt(-1j*omega*sigma*u0*ur);              % Nomalized argument for Bessel function
  J0  =  besselj ( 0 , A.* R );                 % first bessel-function of zeroth order
  J1  =  besselj ( 1 , A.* R );                 % first bessel-function of first order
  z_math     = A./(sigma*2*pi.*R).*J0./J1;      % specific impedance
end