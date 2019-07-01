function nonLinChem(dy,y,p,t)
 dy[1] = -y[1]
 dy[2] = y[1]-(y[2])^2
 dy[3] = (y[2])^2
end
y0 = [1.0;0.0;0.0]
tspan = (0.0,20)
nlc_analytic(u0,p,t) = [exp(-t);
    (2sqrt(exp(-t))besselk(1,2sqrt(exp(-t)))-2besselk(1,2)/besseli(1,2)*sqrt(exp(-t))besseli(1,2sqrt(exp(-t))))/(2besselk(0,2sqrt(exp(-t)))+(2besselk(1,2)/besseli(1,2))besseli(0,2sqrt(exp(-t))));
    -exp(-t)+1+(-2sqrt(exp(-t))*besselk(1,2sqrt(exp(-t)))+sqrt(exp(-t))*besseli(1,2sqrt(exp(-t)))*2besselk(1,2)/besseli(1,2))/(2besselk(0,2sqrt(exp(-t)))+2besselk(1,2)/besseli(1,2)*besseli(0,2sqrt(exp(-t))))]
nonLinChem_f = ODEFunction(nonLinChem,analytic = nlc_analytic)
prob_ode_nonlinchem = ODEProblem(nonLinChem,y0,tspan)
