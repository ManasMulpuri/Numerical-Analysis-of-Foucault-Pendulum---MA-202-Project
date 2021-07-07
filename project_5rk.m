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
x_actual = zeros(time/h+1,0);
y_actual = zeros(time/h+1,0);
m = zeros(time/h+1,0);
n = zeros(time/h+1,0);
t_vector = zeros(time/h+1,0);
x(1) = 1;
y_actual(1) = -0.013604;
x_actual(1) = 1;
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
    
    a2 = dx_dt(x(i)+a1*h/4,y(i)+b1*h/4,m(i)+c1*h/4,n(i)+d1*h/4,omega,phi,wo,r);
    b2 = dy_dt(x(i)+a1*h/4,y(i)+b1*h/4,m(i)+c1*h/4,n(i)+d1*h/4,omega,phi,wo,r);
    c2 = dm_dt(x(i)+a1*h/4,y(i)+b1*h/4,m(i)+c1*h/4,n(i)+d1*h/4,omega,phi,wo,r);
    d2 = dn_dt(x(i)+a1*h/4,y(i)+b1*h/4,m(i)+c1*h/4,n(i)+d1*h/4,omega,phi,wo,r);
    
    a3 = dx_dt(x(i)+a1*h/8+a2*h/8,y(i)+b1*h/8+b2*h/8,m(i)+c1*h/8+c2*h/8,n(i)+d1*h/8+d2*h/8,omega,phi,wo,r);
    b3 = dy_dt(x(i)+a1*h/8+a2*h/8,y(i)+b1*h/8+b2*h/8,m(i)+c1*h/8+c2*h/8,n(i)+d1*h/8+d2*h/8,omega,phi,wo,r);
    c3 = dm_dt(x(i)+a1*h/8+a2*h/8,y(i)+b1*h/8+b2*h/8,m(i)+c1*h/8+c2*h/8,n(i)+d1*h/8+d2*h/8,omega,phi,wo,r);
    d3 = dn_dt(x(i)+a1*h/8+a2*h/8,y(i)+b1*h/8+b2*h/8,m(i)+c1*h/8+c2*h/8,n(i)+d1*h/8+d2*h/8,omega,phi,wo,r);
    
    a4 = dx_dt(x(i)-a2*h/2+a3*h,y(i)-b2*h/2+b3*h,m(i)-c2*h/2+c3*h,n(i)-d2*h/2+d3*h,omega,phi,wo,r);
    b4 = dy_dt(x(i)-a2*h/2+a3*h,y(i)-b2*h/2+b3*h,m(i)-c2*h/2+c3*h,n(i)-d2*h/2+d3*h,omega,phi,wo,r);
    c4 = dm_dt(x(i)-a2*h/2+a3*h,y(i)-b2*h/2+b3*h,m(i)-c2*h/2+c3*h,n(i)-d2*h/2+d3*h,omega,phi,wo,r);
    d4 = dn_dt(x(i)-a2*h/2+a3*h,y(i)-b2*h/2+b3*h,m(i)-c2*h/2+c3*h,n(i)-d2*h/2+d3*h,omega,phi,wo,r);
    
    a5 = dx_dt(x(i)+a1*h*3/16+a4*h*9/16,y(i)+b1*h*3/16+b4*h*9/16,m(i)+c1*h*3/16+c4*h*9/16,n(i)+d1*h*3/16+d4*h*9/16,omega,phi,wo,r);
    b5 = dy_dt(x(i)+a1*h*3/16+a4*h*9/16,y(i)+b1*h*3/16+b4*h*9/16,m(i)+c1*h*3/16+c4*h*9/16,n(i)+d1*h*3/16+d4*h*9/16,omega,phi,wo,r);
    c5 = dm_dt(x(i)+a1*h*3/16+a4*h*9/16,y(i)+b1*h*3/16+b4*h*9/16,m(i)+c1*h*3/16+c4*h*9/16,n(i)+d1*h*3/16+d4*h*9/16,omega,phi,wo,r);
    d5 = dn_dt(x(i)+a1*h*3/16+a4*h*9/16,y(i)+b1*h*3/16+b4*h*9/16,m(i)+c1*h*3/16+c4*h*9/16,n(i)+d1*h*3/16+d4*h*9/16,omega,phi,wo,r);
    
    a6 = dx_dt(x(i)-a1*h*3/7+a2*h*2/7+a3*h*12/7-a4*h*12/7+a5*h*8/7,y(i)-b1*h*3/7+b2*h*2/7+b3*h*12/7-b4*h*12/7+b5*h*8/7,m(i)-c1*h*3/7+c2*h*2/7+c3*h*12/7-c4*h*12/7+c5*h*8/7,n(i)-d1*h*3/7+d2*h*2/7+d3*h*12/7-d4*h*12/7+d5*h*8/7,omega,phi,wo,r);
    b6 = dy_dt(x(i)-a1*h*3/7+a2*h*2/7+a3*h*12/7-a4*h*12/7+a5*h*8/7,y(i)-b1*h*3/7+b2*h*2/7+b3*h*12/7-b4*h*12/7+b5*h*8/7,m(i)-c1*h*3/7+c2*h*2/7+c3*h*12/7-c4*h*12/7+c5*h*8/7,n(i)-d1*h*3/7+d2*h*2/7+d3*h*12/7-d4*h*12/7+d5*h*8/7,omega,phi,wo,r);
    c6 = dm_dt(x(i)-a1*h*3/7+a2*h*2/7+a3*h*12/7-a4*h*12/7+a5*h*8/7,y(i)-b1*h*3/7+b2*h*2/7+b3*h*12/7-b4*h*12/7+b5*h*8/7,m(i)-c1*h*3/7+c2*h*2/7+c3*h*12/7-c4*h*12/7+c5*h*8/7,n(i)-d1*h*3/7+d2*h*2/7+d3*h*12/7-d4*h*12/7+d5*h*8/7,omega,phi,wo,r);
    d6 = dn_dt(x(i)-a1*h*3/7+a2*h*2/7+a3*h*12/7-a4*h*12/7+a5*h*8/7,y(i)-b1*h*3/7+b2*h*2/7+b3*h*12/7-b4*h*12/7+b5*h*8/7,m(i)-c1*h*3/7+c2*h*2/7+c3*h*12/7-c4*h*12/7+c5*h*8/7,n(i)-d1*h*3/7+d2*h*2/7+d3*h*12/7-d4*h*12/7+d5*h*8/7,omega,phi,wo,r);
    
    x(i+1) = x(i) + (1/90)*(7*a1+32*a3+12*a4+32*a5+7*a6)*h;
    y(i+1) = y(i) + (1/90)*(7*b1+32*b3+12*b4+32*b5+7*b6)*h;
    m(i+1) = m(i) + (1/90)*(7*c1+32*c3+12*c4+32*c5+7*c6)*h;
    n(i+1) = n(i) + (1/90)*(7*d1+32*d3+12*d4+32*d5+7*d6)*h;
    
    
    t_vector(i+1) = t_vector(i)+h;
    t = t+ h;
    x_actual(i+1) = 0.5000504*cos(1.1073*t)+0.499949*cos(1.1074*t);
    y_actual(i+1) = 0.5000504*sin(1.1073*t)+0.499949*sin(1.1074*t)-0.013604;
    i = i+1;
end
%ans = abs(x-x_actual);
plot(x,y,'g')
xlabel('x')
ylabel('y')
title('X vs Y')
hold off
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
    d = -2*omega*m*sin(phi) - wo*wo*y - r*omega*omega*cos(phi)*sin(phi);
end