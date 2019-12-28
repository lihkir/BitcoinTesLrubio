%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=convec(u, h, K)
%int:convec_type
%double:h
%matrix:u,K

global convec_type;

if ( convec_type == 0 )
  K=convec_pvm(u, h, K);
else
  K=convec_glf(u, h, K);
end