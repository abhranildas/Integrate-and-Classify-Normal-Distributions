f=@(x) (x-1).*(x-3-1i).*(x-3+1i).*(x-4);
fplot(f);
ylim([-1 100])
refline(0,0)

%%
x1=1;
x2=3;

x3=2;
x4=4;

p1=poly([x1 x2 x3 x4]);
p2=poly([x1 x4]);

%fplot(@(x) (x-x1).*(x-x2));
hold on

fplot(@(x) polyval(p1,x));
fplot(@(x) polyval(p2,x));
fplot(@(x) polyval(p1,x)+polyval(p2,x));
ylim([-5 5])
refline(0,0)


a=(x2+x3)/2; b=(x3-x2)/2*sqrt(sign(x3-x2));
x2=a+b;
x3=a-b;

p3=poly([x1 x2 x3 x4])

fplot(@(x) polyval(p3,x));
ylim([-4 4])
refline(0,0)
hold off

%% check what the different quartic coefficients do
q4=1;
q3=1;
q2=1;
q1=1;
q0=1;
hold on;
for q3=1:4
    f=@(x) q4*x.^4+q3*x.^3+q2*x.^2+q1*x+q0;
    fplot(f,'color','k');
end


%%
fsurf(@(x,y) x.*y,'ShowContours','on')

fsurf(@(x,y) x.^2-y.^2.*sign(y),'ShowContours','on')

fsurf(@(x,y) (x+y).^2-(x-y).^2.*sign(x-y),'ShowContours','on')

%% implicit function for the intersection of two quadratic domains
quad1.q2=[1 1; 1 1];
quad1.q1=[-1;0];
quad1.q0=-1;

quad2.q2=[1 -1; -1 1];
quad2.q1=[1;0];
quad2.q0=-1;

f1=quad2fun(quad1);
fimplicit(f1)
f2=quad2fun(quad2);
f=@(x,y) min(f1(x,y),f2(x,y));

fimplicit(f)
