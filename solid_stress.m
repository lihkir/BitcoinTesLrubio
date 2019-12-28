function v = solid_stress(phi)

global diff_idx;

if (diff_idx == 1)
    v = sige(phi);
elseif (diff_idx == 2)
    v = sige_reg(phi);
else 
    error('Undefined diff_idx!!');
end