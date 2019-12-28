function Sa=minmax_charspeed(ua, Sa)

%int: M2, N2, i, j, gc
%vector: uah, uai, ul, ur, uh
%matrix: ua, Sa

global gc;

M2=size(ua, 2)-2*gc;
N2=size(ua, 1);

ul=zeros(N2,1);
ur=zeros(N2,1);
uh=zeros(3,1);

for i = gc+1:M2+gc
    for j = 1:N2
        ul(j)=ua(j, i);
        ur(j)=ua(j, i+1);
    end
    uh=lminmax_charspeed(ul, ur, uh);
    for j = 1:2
        Sa(j, i)=uh(j);
    end
end