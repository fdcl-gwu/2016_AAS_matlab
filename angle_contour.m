% 5 May 16
% angle between the two gradient vectors


function [] = angle_contour(fig_handle,constants)


plane = constants.plane; % plane to plot ZVC curves
pot_model = constants.pot_model;
num_f = constants.asteroid_grav.num_f;
%Contour plot of potential
m = -2:0.01:2;
[X,Y]=meshgrid(m);
xmin = -2;xmax = 2;ymin=-2; ymax = 2;zmin = -2;zmax=2;

asteroid_grav = constants.asteroid_grav;
% zvc levels

contour_levels = [0:1:10];
switch plane
    case 'xy' % X-Y
        theta = zeros(size(X));
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);Y(ii,jj);0];
                
                switch num_f
                    case 1024
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                    case 4092
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                end
                [~,Ugm] = mascon_potential(state, asteroid_grav,constants);
                
                angle = asin(norm(cross(Ugm,Ugp))/norm(Ugp)/norm(Ugm));
                theta(ii,jj) = angle;
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        
        hold all
        grid on
        axis([xmin xmax ymin ymax zmin zmax])
        
        %         title_string = ['4769 Castalia Energy Contour'];
        %         title(title_string,'interpreter','latex','FontName','times','FontSize',18)
        xlabel('X Axis ','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Y Axis ','interpreter','latex','FontName','times','FontSize',18)
        %                 zlabel('Z Axis ','interpreter','latex','FontName','times','FontSize',18)
        %
        
        contour(X,Y,theta*180/pi,contour_levels,'r', 'ShowText','on');
        
        
        
        
    case 'xz' %x-Z
        U = zeros(size(X));
        U_grad = zeros(size(X));
        
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);0;Y(ii,jj)];
                switch num_f
                    case 1024
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                    case 4092
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                end
                [~,Ugm] = mascon_potential(state, asteroid_grav,constants);
                
                angle = asin(norm(cross(Ugm,Ugp))/norm(Ugp)/norm(Ugm));
                theta(ii,jj) = angle;
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        
        axis([xmin xmax ymin ymax zmin zmax])
        xlabel('X Axis ','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z Axis ','interpreter','latex','FontName','times','FontSize',18)
         contour(X,Y,theta*180/pi,contour_levels,'r', 'ShowText','on');
        
        
    case 'yz' % Y-Z
        U = zeros(size(X));
        U_grad = zeros(size(X));
        
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [0;X(ii,jj);Y(ii,jj)];
                switch num_f
                    case 1024
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                    case 4092
                        [~,Ugp, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                end
                [~,Ugm] = mascon_potential(state, asteroid_grav,constants);
                
                angle = asin(norm(cross(Ugm,Ugp))/norm(Ugp)/norm(Ugm));
                theta(ii,jj) = angle;
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        xlabel('Y Axis ','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z Axis ','interpreter','latex','FontName','times','FontSize',18)
        axis([xmin xmax ymin ymax zmin zmax])
         contour(X,Y,theta*180/pi,contour_levels,'r', 'ShowText','on');
        
        
        
end

end