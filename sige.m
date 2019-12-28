function v = sige(phi)

global phic sigma0 kexp;

if (phi <= phic)
    v = 0;
else
    v = sigma0*((phi/phic)^kexp - 1);
end