function v = exp_reg(phi)

global epsilon phic;

v = exp(-epsilon/(phi - phic)^2);