
function rotation_matrix = PQW2ECI(arg_prg, inc_angle, RAAN)

% arg_prg = deg2rad(arg_prg);
% inc_angle = deg2rad(inc_angle);
% RAAN = deg2rad(RAAN);

dcm_w3 = [cos(arg_prg) sin(arg_prg) 0 ; -sin(arg_prg) cos(arg_prg) 0 ; 0 0 1] ;         %R(w,3'')
dcm_i1 = [1 0 0 ; 0 cos(inc_angle) sin(inc_angle) ; 0 -sin(inc_angle) cos(inc_angle)] ; %R(i,1')
dcm_ohm3 = [cos(RAAN) sin(RAAN) 0 ; -sin(RAAN) cos(RAAN) 0 ; 0 0 1] ;                   %R(ohm,3)

rotation_matrix = (dcm_w3 * dcm_i1 * dcm_ohm3)';

end
