% 4 May 16
% Zero Velocity Curves for a range of jacobi energy values

function [] = zvc_contour(fig_handle,constants)

plane = constants.plane; % plane to plot ZVC curves
pot_model = constants.pot_model;
num_f = constants.asteroid_grav.num_f;

%Contour plot of potential
[X,Y]=meshgrid(-2:0.01:2);
xmin = -2;xmax = 2;ymin=-2; ymax = 2;zmin = -2;zmax=2;

asteroid_grav = constants.asteroid_grav;
% zvc levels J should be a vector of lines where contour should be plotted
scale = 1e6; % units are in km^2/sec^2
zvc_lvls = 0.01

switch plane
    case {'xy','XY'} % X-Y
        U = zeros(size(X));
        U_grad = zeros(size(X));
        
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);Y(ii,jj);0];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U(ii,jj),Ug] = mascon_potential(state, asteroid_grav,constants);
                end
                U_grad(ii,jj) = norm(Ug);
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        
        hold all
        grid on
        axis([xmin xmax ymin ymax zmin zmax])
        
        % contour(X,Y,U_grad*1e6,[0.05:0.05:0.25],'r', 'ShowText','on'); % grav attraction
        V = 1/2.*constants.omega^2.*(X.^2+Y.^2) + U;
        contour(X,Y,V*scale,zvc_lvls,'r', 'ShowText','on');
        
        xlabel('X (km)','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Y (km)','interpreter','latex','FontName','times','FontSize',18)
        
    case {'xz','XZ'} %x-Z
        U = zeros(size(X));
        U_grad = zeros(size(X));
        Ulaplace = zeros(size(X));
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);0;Y(ii,jj)];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U(ii,jj),Ug] = mascon_potential(state, asteroid_grav,constants);
                end
                U_grad(ii,jj) = norm(Ug);
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        
        axis([xmin xmax ymin ymax zmin zmax])
        
        
        % contour(X,Y,U_grad*1e6,[0.05:0.05:0.25],'r', 'ShowText','on');
        V = 1/2.*constants.omega^2.*(X.^2+zeros(size(X)).^2) + U;
        contour(X,Y,V*scale,zvc_lvls,'r', 'ShowText','on');
        
        xlabel('X (km)','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z (km)','interpreter','latex','FontName','times','FontSize',18)
    case {'yz','YZ'} % Y-Z
        U = zeros(size(X));
        U_grad = zeros(size(X));
        Ulaplace = zeros(size(X));
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [0;X(ii,jj);Y(ii,jj)];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U(ii,jj),Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U(ii,jj),Ug] = mascon_potential(state, asteroid_grav,constants);
                end
                U_grad(ii,jj) = norm(Ug);
                
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        
        axis([xmin xmax ymin ymax zmin zmax])
        
        
        % contour(X,U,U_grad*1e6,[0.05:0.05:0.25],'r', 'ShowText','on');
        V = 1/2.*constants.omega^2.*(zeros(size(Y)).^2+X.^2) + U;
        contour(X,Y,V*scale,zvc_lvls,'r', 'ShowText','on');
        
        xlabel('Y (km)','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z (km)','interpreter','latex','FontName','times','FontSize',18)
end

end