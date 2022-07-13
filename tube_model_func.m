% R_DC tubemodell (function)
function  R_dc = tube_model_func( delta ,r_outer, sigma )
  R_rod   = pi/2/sigma*r_outer^2;                       % rod resistance
  r_inner = max(r_outer - delta ,0);                    % inner radius of tube belong to skindepth
  R_tube  = pi./2./sigma.*(r_outer.^2 - r_inner.^2) ;   % tube resistance
  R_dc    = R_rod./R_tube;                              % Calculated DC-Resistance
end