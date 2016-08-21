% 3 May 16
% mascon gravity model
function [U,Ugrad,U_grad_mat] = mascon_potential(state,asteroid_grav,constants)

pos = reshape(state(1:3),1,3);
M = constants.M;
G = constants.G;
% get all the vectors for the faces
V = asteroid_grav.V;
num_v = asteroid_grav.num_v;

% calculate the mass of each particle
m = M/num_v;
% compute the distance from the state to each particle
r = repmat(pos,num_v,1)-V;

rnorm = sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2);
rnorm3 = rnorm.^3;
rnorm5 = rnorm.^5;


% sum up potential and gradient
U = G*sum(m./rnorm);
Ugrad = -G*sum(m.*r./repmat(rnorm3,1,3))';

% calculate second order partial matrix
Uxx = sum( m.*(-3.*r(:,1).^2./rnorm5 + 1./rnorm3));
Uyy = sum( m.*(-3.*r(:,2).^2./rnorm5 + 1./rnorm3));
Uzz = sum( m.*(-3.*r(:,3).^2./rnorm5 + 1./rnorm3));

Uxy = sum( m.*(-3.*r(:,1).*r(:,2)./rnorm5 )); Uyx = Uxy;
Uzx = sum( m.*(-3.*r(:,3).*r(:,1)./rnorm5 )); Uxz = Uzx;
Uyz =sum( m.*(-3.*r(:,3).*r(:,2)./rnorm5 )); Uzy = Uyz;

U_grad_mat = -G*[Uxx Uxy Uxz;...
             Uyx Uyy Uyz;...
             Uzx Uzy Uzz];