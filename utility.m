function [ct, at, ht, ut] = utility(par, gr, b, theta, w_l, w_h, r, P_h)

a_opt = @(a, h) (b - a + w_l*(1-h) - P_h*h)^(-par.psi) - ...
                par.beta*(1+r)*( ((h*theta)^par.eta)    *(w_h + (1+r)*a)^(-par.psi) + ...
                             (1-((h*theta)^par.eta))*(w_l + (1+r)*a)^(-par.psi));

% lim_min = -w_l/(1+r);
% lim_max_0 = b + w_l ;
% lim_max_1 = b - P_h ;
% 
% a_0 = bisection(@(x) a_opt(x,0), lim_min, lim_max_0);
% a_1 = bisection(@(x) a_opt(x,1), lim_min, lim_max_1);
                         
a_0 = fsolve(@(x) a_opt(x,0), 0);   % Without college education
a_1 = fsolve(@(x) a_opt(x,1), 0);   % With college education

% Savings must be above borrowing constraint, and cannot make consumption
% negative
if a_0 < - gr.Abar
    a_0 = -gr.Abar;
elseif a_0 > b + w_l
    a_0 = b + w_l - eps;
end

if a_1 < - gr.Abar
    a_1 = -gr.Abar;
elseif a_1 > b - P_h
    a_1 = b - P_h - eps;
end

u_0 = ((b - a_0 + w_l)^(1-par.psi))/(1-par.psi) + par.beta * ((w_l + (1+r)*a_0)^(1-par.psi))/(1-par.psi);
u_1 = ((b - a_1 - P_h)^(1-par.psi))/(1-par.psi) + ...
        par.beta*( (theta^par.eta)*((w_h + (1+r)*a_1)^(1-par.psi))/(1-par.psi) + ...
             (1-(theta^par.eta))*((w_l + (1+r)*a_1)^(1-par.psi))/(1-par.psi));

if u_0 >= u_1
    ht = 0;
    at = a_0;
    ut = u_0;
    ct = b - a_0 + w_l;
else
    ht = 1;
    at = a_1;
    ut = u_1;
    ct = b - a_1 - P_h;
end

end