clear all;

omega = 7.292E-5;
phi = 50;
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
T = zeros(time/h+1,0);

x(1) = 1;
y(1) = 0;
m(1) = 1;
n(1) = 0;

i=1;
T(i) = t;

while t<time
    a1 = dx_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    b1 = dy_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    c1 = dm_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    d1 = dn_dt(x(i),y(i),m(i),n(i),omega,phi,wo,r);
    
    a2 = dx_dt(x(i)+a1*h/2,y(i)+b1*h/2,m(i)+c1*h/2,n(i)+d1*h/2,omega,phi,wo,r);
    b2 = dy_dt(x(i)+a1*h/2,y(i)+b1*h/2,m(i)+c1*h/2,n(i)+d1*h/2,omega,phi,wo,r);
    c2 = dm_dt(x(i)+a1*h/2,y(i)+b1*h/2,m(i)+c1*h/2,n(i)+d1*h/2,omega,phi,wo,r);
    d2 = dn_dt(x(i)+a1*h/2,y(i)+b1*h/2,m(i)+c1*h/2,n(i)+d1*h/2,omega,phi,wo,r);
    
    a3 = dx_dt(x(i)+a2*h/2,y(i)+b2*h/2,m(i)+c2*h/2,n(i)+d2*h/2,omega,phi,wo,r);
    b3 = dy_dt(x(i)+a2*h/2,y(i)+b2*h/2,m(i)+c2*h/2,n(i)+d2*h/2,omega,phi,wo,r);
    c3 = dm_dt(x(i)+a2*h/2,y(i)+b2*h/2,m(i)+c2*h/2,n(i)+d2*h/2,omega,phi,wo,r);
    d3 = dn_dt(x(i)+a2*h/2,y(i)+b2*h/2,m(i)+c2*h/2,n(i)+d2*h/2,omega,phi,wo,r);
    
    a4 = dx_dt(x(i)+a3*h,y(i)+b3*h,m(i)+c3*h,n(i)+d3*h,omega,phi,wo,r);
    b4 = dy_dt(x(i)+a3*h,y(i)+b3*h,m(i)+c3*h,n(i)+d3*h,omega,phi,wo,r);
    c4 = dm_dt(x(i)+a3*h,y(i)+b3*h,m(i)+c3*h,n(i)+d3*h,omega,phi,wo,r);
    d4 = dn_dt(x(i)+a3*h,y(i)+b3*h,m(i)+c3*h,n(i)+d3*h,omega,phi,wo,r);
    
   
    x(i+1) = x(i) + (1/6)*(a1+2*a2+2*a3+a4)*h;
    y(i+1) = y(i) + (1/6)*(b1+2*b2+2*b3+b4)*h;
    m(i+1) = m(i) + (1/6)*(c1+2*c2+2*c3+c4)*h;
    n(i+1) = n(i) + (1/6)*(d1+2*d2+2*d3+d4)*h;
    
    t = t+ h;
    i = i+1;
    T(i) = t;
end

plot(T,x)
xlabel('t')
ylabel('x')
title('Movement of Pendulum in x vs t')


function a = dx_dt(x,y,m,n,omega,phi,wo,r)
    a = m;
end

function b = dy_dt(x,y,m,n,omega,phi,wo,r)
    b = n;
end

function c = dm_dt(x,y,m,n,omega,phi,wo,r)
    c = 2*omega*n*sind(phi) - wo*wo*x;
end

function d = dn_dt(x,y,m,n,omega,phi,wo,r)
    d = -2*omega*m*sind(phi) - wo*wo*y - r*(omega.^2)*cosd(phi)*sind(phi);
end

