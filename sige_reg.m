function v = sige_reg(phi)

global phic;

if (phi <= phic)
    v = 0;
else
    v = sige(phi)*exp_reg(phi);
end