function y = fullfitfunc_1d(x,mu_eff, D, a)
l = D*3;
z0 = l;
z1 = -z0 - 4*D;
x0 = abs(x-z0);
x1 = abs(x-z1);
y = a*(exp(-mu_eff.*x0) - exp(-mu_eff.*x1));
end
