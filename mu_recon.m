function mu = mu_recon(y,x,mu_eff, D)
l = D*3;
z0 = l;
z1 = -z0 - 4*D;
first_coef = exp(mu_eff*z0) - exp(mu_eff*z0);
second_coef = z1*exp(mu_eff*z0) - z0*exp(mu_eff*z0);
mu = 4*pi*D*(x.^2 -(z0+z1).*x +z0*z1).*y;
mu = mu./(first_coef.*x - second_coef);
mu = -log(mu)./x;
end
