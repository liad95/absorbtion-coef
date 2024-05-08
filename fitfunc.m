function y = fitfunc(x,a,b,z0,z1)
y = a*(exp(-b*(x-z0)) - exp(-b*(x-z1)))
end