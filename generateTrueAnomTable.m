function [table] = generateTrueAnomTable(n,e,num_nu_pts)

P = 2*pi/n;
t = 0:P/(num_nu_pts-1):P;
M = n.*t;
E = zeros(1,numel(t));
for i = 1:numel(t)
    E(i) = fzero(@(E) E-e.*sin(E)-M(i), M(i),optimoptions('fsolve','Display','off','TolFun',1e-50,'TolX',1e-50));
end%for
nu = acos((cos(E)-e)./(1-e.*cos(E)));
nu(end/2:end) = 2*pi - nu(end/2:end); % True anomaly should be monotonically increasing.

table = [t' nu'];










end%function