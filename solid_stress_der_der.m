function v = solid_stress_der_der(phi)

global diff_idx;

if (diff_idx == 1)
    v = sige_der_der(phi);
elseif (diff_idx == 2)
    v = sige_reg_der_der(phi);
else
    error('Undefined diff_idx!!');
end