function v = sige_der(phi)

global phic sigedc kexp;

if (phi <= phic)
    v = 0;
else
    v = sigedc*((phi/phic)^(kexp - 1));
end