function fh= fhroe(ul, ur, fl, fr)

fl=fn_flux(ul,fl);
fr=fn_flux(ur,fr);

fh=0.5*(fl+fr-abs(fr-fl)*sign(ur-ul));
