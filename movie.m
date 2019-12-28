x=[-2*pi:0.1:2*pi];
fig = figure;
mov = avifile('ejemplo.avi'); % Abrimos el video
for k=1:17 % Para cada uno de los frames
y=sin(x+k*pi/8);
plot(x,y);
F = getframe(fig); % Almacenamos la figura como frame
mov = addframe(mov, F);
end
close(fig);
mov = close(mov); % Cerramos la película

,'compression','None'