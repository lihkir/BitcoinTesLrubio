function v = sige_der_der(phi)

global phic sigeddc kexp;

if (phi <= phic)
    v = 0;
else
    v = sigeddc*((phi/phic)^(kexp - 2));
end