function draw_ellipse(mu,sigma)
t = linspace(0,2*pi) ;
a = 30 ; b = 15 ;
x = a*cos(t) ;
y = b*sin(t) ;
plot(x,y,'r')
axis equal