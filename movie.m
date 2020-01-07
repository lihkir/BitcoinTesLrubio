n=(1:50)*2*pi ;

Y = sin(n*50) ;
hp = plot(Y) ;              %// Generate the initial plot (and retrieve the handle of the graphic object)
ylim([-1,1]) ;              %// Set the Y axes limits (once and for all)

writerObj = VideoWriter('test2.avi'); %// initialize the VideoWriter object
open(writerObj) ;
for t = 1:1000
   Y = sin(n*50/t) ;        %// calculate new Y values
   set(hp,'YData',Y) ;      %// update the plot data (this does not generate a "new" plot), nor resize the axes

   F = getframe ;           %// Capture the frame
   writeVideo(writerObj,F)  %// add the frame to the movie
end
close(writerObj);