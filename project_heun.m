omega = 7.2921150E-5;
phi = 0.872665;
g = 9.81;
l = 8;
wo = sqrt(g/l); 
r = 6371397; %m

time = 2*60*60; %2 hours in sec
h = 0.1;

t = 0;

x = zeros(time/h+1,0);
y = zeros(time/h+1,0);
m = zeros(time/h+1,0);
n = zeros(time/h+1,0);
t_vector = zeros(time/h+1,0);
x(1) = 1;
y(1) = 0;
m(1) = 0;
n(1) = 0;
t_vector(1) = 0;
i=1;

while t<time
    a1 = dx_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    b1 = dy_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    c1 = dm_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    d1 = dn_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    
    a2 = dx_dt(x(i)+a1*h,y(i)+b1*h,m(i)+c1*h,n(i)+d1*h,omega,phi,wo,r);
    b2 = dy_dt(x(i)+a1*h,y(i)+b1*h,m(i)+c1*h,n(i)+d1*h,omega,phi,wo,r);
    c2 = dm_dt(x(i)+a1*h,y(i)+b1*h,m(i)+c1*h,n(i)+d1*h,omega,phi,wo,r);
    d2 = dn_dt(x(i)+a1*h,y(i)+b1*h,m(i)+c1*h,n(i)+d1*h,omega,phi,wo,r);
    
    x(i+1) = x(i) + (a1+a2)*h/2;
    y(i+1) = y(i) + (b1+b2)*h/2;
    m(i+1) = m(i) + (c1+c2)*h/2;
    n(i+1) = n(i) + (d1+d2)*h/2;
    t_vector(i+1) = t_vector(i)+h;
    t = t+ h;
    i = i+1;
end

plot(t_vector,x,'r')
xlabel('t')
ylabel('x')
title('Movement of Pendulum in T vs X')
function a = dx_dt(x,y,m,n,omega,phi,wo,r)
    a = m;
end

function b = dy_dt(x,y,m,n,omega,phi,wo,r)
    b = n;
end

function c = dm_dt(x,y,m,n,omega,phi,wo,r)
    c = 2*omega*n*sin(phi) - wo*wo*x;
end

function d = dn_dt(x,y,m,n,omega,phi,wo,r)
    d = -2*omega*m*sin(phi) - wo*wo*y- r*omega*omega*cos(phi)*sin(phi);
end