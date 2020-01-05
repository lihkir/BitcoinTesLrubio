function v = sige_reg_der(phi)

v = exp_reg(phi)*(sige_der(phi) + 2*sige(phi)*quot_reg(phi));

end