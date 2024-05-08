function f = exp_3d_func(x,y,z,mu_eff, a)
f = a*(exp(-mu_eff.*sqrt(x.^2 + y.^2 + z.^2)));
end
