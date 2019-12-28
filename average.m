function ua=average(ul, ur, s, ua)

%int: i, N
%vector: ul, ur, ua 
%double: s

N = size(ul,2);

for i=1:N
  ua(i)=ul(i)+s*(ur(i)-ul(i));
end
