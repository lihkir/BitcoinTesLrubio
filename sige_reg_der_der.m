function v = sige_reg_der_der(phi)
    
v = exp_reg(phi)*(sige_der_der(phi) + 4*sige_der(phi)*quot_reg(phi) + 2*epsilon*sige(phi)*(2*epsilon - 3*(phi - phic)^2));

end