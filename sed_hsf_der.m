function v=sed_hsf_der(phi)

%double: phi, v, vinf, nexp, phimax, phi0, vphi0, v1phi0

global vinf phimax phi0 nexp;

if (phi < 0 || phi > phimax)
    v = 0;
elseif (phi > phi0)
    vphi0 = (1 - phi0)^(nexp - 2);
    v1phi0 = -(nexp-2)*(1 - phi0)^(nexp - 3);
    v = vinf*(-vphi0 + v1phi0*(1 + phi0) - 2*phi*v1phi0);
else
    v = -(nexp - 1)*vinf*(1 - phi)^(nexp - 2);
end