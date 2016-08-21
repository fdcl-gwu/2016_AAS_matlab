% 19 April 2016
% contour of gravitational attraction or zero velocity curves about 4769 Castalia


function [] = pot_contour(fig_handle,constants)


plane = constants.plane; % plane to plot ZVC curves
pot_model = constants.pot_model;
num_f = constants.asteroid_grav.num_f;
%Contour plot of potential
m = -2:0.01:2;
[X,Y]=meshgrid(m);
xmin = -2;xmax = 2;ymin=-2; ymax = 2;zmin = -2;zmax=2;

asteroid_grav = constants.asteroid_grav;
% zvc levels
scale = 1e6;
contour_levels = [0.001:0.005:0.25];
switch plane
    case {'xy','XY'} % X-Y
        Ugrid = zeros(size(X));
        U_grad = zeros(size(X));

        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);Y(ii,jj);0];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U,Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U,Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U,Ug] = mascon_potential(state, asteroid_grav,constants);
                    case 'diff'
                        % calculate norm of difference of the attraction
                        % models
                        switch num_f
                            case 1024
                                [Up,~, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [Up,~, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                        [Um,~] = mascon_potential(state, asteroid_grav,constants);
                        U = abs(Up-Um);        
                        scale = 1e7;
                        
                        
                end
                
                Ugrid(ii,jj) = U;
               
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        
        hold all
        grid on
        axis([xmin xmax ymin ymax zmin zmax])
        
        %         title_string = ['4769 Castalia Energy Contour'];
        %         title(title_string,'interpreter','latex','FontName','times','FontSize',18)
                xlabel('X (km)','interpreter','latex','FontName','times','FontSize',18)
                ylabel('Y (km)','interpreter','latex','FontName','times','FontSize',18)
%                 zlabel('Z Axis ','interpreter','latex','FontName','times','FontSize',18)
        %
        
        contour(X,Y,Ugrid*scale,contour_levels,'r', 'ShowText','on'); % grav attraction
    
        
        
        
    case {'xz','XZ'} %x-Z
        Ugrid = zeros(size(X));
        U_grad = zeros(size(X));
        
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [X(ii,jj);0;Y(ii,jj)];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U,Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U,Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U,Ug] = mascon_potential(state, asteroid_grav,constants);
                    case 'diff'
                        
                        % calculate norm of difference of the attraction
                        % models
                        switch num_f
                            case 1024
                                [Up,~, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [Up,~, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                        [Um,~] = mascon_potential(state, asteroid_grav,constants);
                        U = abs(Up-Um);      
                        scale = 1e7;
                        
                end
                
                Ugrid(ii,jj) = U;
               
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        
        axis([xmin xmax ymin ymax zmin zmax])
        xlabel('X (km)','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z (km)','interpreter','latex','FontName','times','FontSize',18)
        contour(X,Y,Ugrid*scale,contour_levels,'r', 'ShowText','on'); % grav attraction
        
        
    case {'yz','YZ'} % Y-Z
        Ugrid = zeros(size(X));
        U_grad = zeros(size(X));
        
        
        % loop over the grid
        for ii = 1:length(X)
            for jj = 1:length(Y)
                state = [0;X(ii,jj);Y(ii,jj)];
                switch pot_model
                    case 'polyhedron'
                        switch num_f
                            case 1024
                                [U,Ug, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [U,Ug, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                    case 'mascon'
                        [U,Ug] = mascon_potential(state, asteroid_grav,constants);
                    case 'diff'
                        
                        % calculate norm of difference of the attraction
                        % models
                        switch num_f
                            case 1024
                                [Up,~, ~, ~] = polyhedron_potential_mex_1024(state, asteroid_grav);
                            case 4092
                                [Up,~, ~, ~] = polyhedron_potential_mex_4092(state, asteroid_grav);
                        end
                        [Um,~] = mascon_potential(state, asteroid_grav,constants);
                        U = abs(Up-Um);      
                        scale = 1e7;
                        
                end
                
                Ugrid(ii,jj) = U;
               
            end
        end
        
        % plot contours about the lagrange points
        set(0,'CurrentFigure',fig_handle)
        hold all
        grid on
        xlabel('Y (km)','interpreter','latex','FontName','times','FontSize',18)
        ylabel('Z (km)','interpreter','latex','FontName','times','FontSize',18)
        axis([xmin xmax ymin ymax zmin zmax])
        contour(X,Y,Ugrid*scale,contour_levels,'r', 'ShowText','on'); % grav attraction
        
        
        
end

end