function m=minmod(a, b)

%double: m, a, b

if ( a*b < 0 )
  m=0;
else
  m=min(abs(a), abs(b));
  if ( a < 0 )
    m=-m;
  end
end