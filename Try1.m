%% Try to plot  the residual function a_opt

%1. Defining the handle function to solve

%For the parameters, run the file "initialization.m".
%Rest of parameters not defined in the file
theta=0.4;
b=2;
h=1;
Ph=1;

a_opt = @(a, h) (b - a + w_l*(1-h) - P_h*h)^(-par.psi) - ...
                par.beta*(1+r)*( ((h*theta)^par.eta)    *(w_h + (1+r)*a)^(-par.psi) + ...
                             (1-((h*theta)^par.eta))*(w_l + (1+r)*a)^(-par.psi));
                         


%Plot the function for a grid of a from min to max
amin=-w_l/(1+r)+1.0e-2;
amax=b+w_l*(1-h)-Ph*h-1.0e-2;
agrid=amin:0.001:amax;
[~,MM]=size(agrid);
yyyy=zeros(MM,1);
for ii=1:1:MM
    yyyy(ii)=a_opt(agrid(ii),1);
end
plot(agrid,yyyy)

wgrid = 3.0:0.002:3.05;
rgrid = 0.04:0.001:0.05;

mat = zeros(numel(wgrid),numel(rgrid));
for iw_h = 1:numel(wgrid)
    for ir = 1:numel(rgrid)
        mat(iw_h, ir) = general_eq(par, gr, P_h, w_l, wgrid(iw_h), rgrid(ir));
    end
end








% This graphs the region of the people who decide to study

mat = zeros(gr.nb, gr.ntheta);
tic
for ib = 1:gr.nb
    for itheta = 1:gr.ntheta
        [ct, at, ht, ut] = utility(par, gr, gr.bgrid(ib), gr.thetagrid(itheta), w_l, w_h, r, P_h);
        mat(ib, itheta) = ht;
    end
end
toc
v=[1,1];

figure(2);
contourf(gr.thetagrid, gr.bgrid, mat,v)
xlabel('\theta'); ylabel('Bequests')
title(['People who study'])

% 
