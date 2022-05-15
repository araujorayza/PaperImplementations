x0 = [0.02;0.1]
p = 1;
tspan=[0,5];
[t,x] = ode45(@(t,x) fp(t,x,p),tspan,x0);
plot(t,x) %one of the states grows unbounded when p=1









