function[result]=capillary_wave_complex_analytical_solution_Prosperetti_1981(s,sigma,nu_u,nu_l,rho_u,rho_l,k_x,g,a_0,u_0,k_z)
k_z=0;
k=(k_x^2+k_z^2)^0.5;
omega_0=(((rho_l-rho_u)/(rho_l+rho_u))*g*k+(sigma/(rho_l+rho_u))*k^3)^0.5;
mu_u=nu_u*rho_u;
mu_l=nu_l*rho_l;
lambda_u=(ones(size(s))*(k^2)+s/nu_u).^0.5;
lambda_l=(ones(size(s))*(k^2)+s/nu_l).^0.5;

Lambda_i_j=4*k*(-rho_l*rho_u*s+k*(mu_u-mu_l)*(rho_u*(ones(size(s))*k-lambda_l)-rho_l*(ones(size(s))*k-lambda_u))+k^2*(mu_l-mu_u)^2*(ones(size(s))*k-lambda_l).*(ones(size(s))*k-lambda_u).*s.^(-1))./((rho_l+rho_u)*(rho_l*(ones(size(s))*k-lambda_u)+rho_u*(ones(size(s))*k-lambda_l)));
result=(1./s).*(ones(size(s))*a_0+(s*u_0-ones(size(s))*omega_0^2*a_0)./(s.^2+Lambda_i_j.*s+ones(size(s))*omega_0^2));