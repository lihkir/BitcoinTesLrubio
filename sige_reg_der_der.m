function v = sige_reg_der_der(phi)

global phic epsilon;

if (phi <= phic)
    v = 0;
else
    v = sige_der_der(phi)*exp_reg(phi) + 4*sige_der(phi)*exp_reg(phi)*quot_reg(phi) + 2*epsilon*sige(phi)*exp_reg(phi)*(2*epsilon - 3*(phi - phic)^2);
end