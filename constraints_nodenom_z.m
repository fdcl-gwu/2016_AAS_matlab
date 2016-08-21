function [m, dmdx] = constraints_nodenom_z(statef,phi_d,constants)
% terminal constraints that should be satisfied
num_con = constants.num_con;
xf = statef(1);
yf = statef(2);
zf = statef(3);

xdf = statef(4);
ydf = statef(5);
zdf = statef(6);

xcf = constants.xcf;
xcen = constants.center_vec;

pmap = constants.pmap; % vector of poincare map indicies

% form difference states
x1 = xf - xcf(pmap(1));
x2 = zf - xcf(pmap(2));
x3 = xdf - xcf(pmap(3));
x4 = zdf - xcf(pmap(4));

% extract out the angles
phi1 = phi_d(1);
phi2 = phi_d(2);
phi3 = phi_d(3);

m = zeros(num_con,1);
dmdx = zeros(4,6);

m(1) = yf;
m(2) = zf;
m(3) = sin(phi1)^2*(x1^2+x2^2+x3^2)-x1^2;
m(4) = sin(phi2)^2*(x2^2+x3^2+x4^2)- x2^2;
m(5) = sin(phi3)^2*(2*x3^2+2*x3*sqrt(x4^2+x3^2)+2*x4^2)-x3-sqrt(x4^2+x3^2);

dxidx = eye(6,6);
dxidx = dxidx(pmap,:);

dmdxi = [0,0,0,0;...
    0,1,0,0;...
    2*sin(phi1)^2*[x1-2*x1,x2,x3,x4];...
    2*sin(phi2)^2*[0,x2-2*x2,x3,x4];...
    [0,0,sin(phi3)^2*(4*x3+2*sqrt(x4^2+x3^2)+2*x3^2*sqrt(x4^2+x3^2))-1-x3*sqrt(x4^2+x3^2),sin(phi3)^2*(2*x3*x4*sqrt(x4^2+x3^2)+4*x4)-x4*sqrt(x4^2+x3^2)]];

dmdx = dmdxi*dxidx;

end