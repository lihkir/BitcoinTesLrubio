function v = sige_reg_der(phi)

global phic;

if (phi <= phic)
    v = 0;
else
    v = sige_der(phi)*exp_reg(phi) + 2*sige(phi)*exp_reg(phi)*quot_reg(phi);
end