% 19 July 2016
% input angles for 4-D Poincare section

function [phi_d] = poincare_angles(index_vec,num_angles)

ii = index_vec(1);
jj = index_vec(2);
kk = index_vec(3);

phi1 = linspace(0,180,num_angles)'*pi/180;
phi2 = phi1;
phi3 = linspace(0,360,num_angles)'*pi/180;

phi_d = [phi1(ii) phi2(jj) phi3(kk)];