function v=sed_hsf(phi)

%double: phi, v, vinf, nexp, phimax, phi0, vphi0, v1phi0

global vinf phimax phi0 nexp;

if (phi < 0 || phi > phimax)
    v = 0;
elseif (phi > phi0)
    vphi0 = (1 - phi0)^(nexp - 2);
    v1phi0 = -(nexp - 2)*(1 - phi0)^(nexp - 3);
    v = vinf*(1 - phi)*(vphi0 + v1phi0*(phi - phi0)); 
else
    v = vinf*(1 - phi)^(nexp - 1);
end
