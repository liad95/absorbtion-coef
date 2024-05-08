function y = singleexpfitfunc(x,mu_eff, l, a)
x0 = abs(x-l);
y = a*exp(-mu_eff.*x0)./x0;
end
