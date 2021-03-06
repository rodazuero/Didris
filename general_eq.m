function error = general_eq(par, gr, P_h, w_l, w_h, r)

H = zeros(gr.nb, gr.ntheta);
L = zeros(gr.nb, gr.ntheta);
A = zeros(gr.nb, gr.ntheta);

for ib = 1:gr.nb
    for itheta = 1:gr.ntheta
        [ct, at, ht, ut] = utility(par, gr, gr.bgrid(ib), gr.thetagrid(itheta), w_l, w_h, r, P_h);
        H(ib, itheta) = ht;
        L(ib, itheta) = 1 - ht;
        A(ib, itheta) = at;
    end
end

HH = sum(sum(H));
LL = sum(sum(L));
KK = sum(sum(A));

w_l_res = par.A * par.gamma * (par.alpha_h * par.z * (HH^(1-(1/par.phi))) + (1 - par.alpha_h) * (LL^ (1-(1/par.phi))))^(((par.gamma * par.phi)/(par.phi - 1))-1) ...
            * (1-par.alpha_h)*(LL^(-1/par.phi)) * (KK^(1-par.gamma));
w_h_res = par.A * par.gamma * (par.alpha_h * par.z * (HH^(1-(1/par.phi))) + (1 - par.alpha_h) * (LL^ (1-(1/par.phi))))^(((par.gamma * par.phi)/(par.phi - 1))-1) ...
            * par.alpha_h * par.z * (HH^(-1/par.phi)) * (KK^(1-par.gamma));
r_res   = par.A * (par.alpha_h * par.z * (HH^(1-(1/par.phi))) + (1 - par.alpha_h) * (LL^ (1-(1/par.phi))))^((par.gamma * par.phi)/(par.phi - 1)) ...
            * (1-par.gamma) * (KK^(-par.gamma));

error = (w_l - w_l_res)^2 + (w_h - w_h_res)^2 + (r - r_res)^2;
        
end