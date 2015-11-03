

w_l=1; w_h = ; r=0.05;

error = general_eq(par, gr, P_h, w_l, w_h, r)





equi = fminunc(@(x) general_eq(par, gr, P_h, x(1), x(2), x(3)), [w_l, w_h, r])

wgrid = 3.0:0.002:3.05;
rgrid = 0.04:0.005:0.05;

mat = zeros(numel(wgrid),numel(rgrid));
for iw_h = 1:numel(wgrid)
    for ir = 1:numel(rgrid)
        mat(iw_h, ir) = general_eq(par, gr, P_h, w_l, wgrid(iw_h), rgrid(ir));
    end
end








% This graphs the region of the people who decide to study

mat = zeros(gr.nb, gr.ntheta);

for ib = 1:gr.nb
    for itheta = 1:gr.ntheta
        [ct, at, ht, ut] = utility(par, gr, gr.bgrid(ib), gr.thetagrid(itheta), w_l, w_h, r, P_h);
        mat(ib, itheta) = ht;
    end
end

v=[1,1];

figure(2);
contourf(gr.thetagrid, gr.bgrid, mat,v)
xlabel('\theta'); ylabel('Bequests')
title(['People who study'])

% 