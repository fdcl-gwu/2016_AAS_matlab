% events cross function
function [value,isterminal,direction]=events_xcross(t,x,constants)

% logic to end integration based on type of orbit desired
direction=0; % any direction is valid
switch constants.periodic_diffcorr_section
    case 'x_axis'
        
        if abs(t) > 100 && x(1) > 0 % wait a short time before checking to avoid initial crossing
            isterminal = 1;
        else
            isterminal = 0;
        end
        value=x(2);  % stop at x axis cross (y is 0 or x(2) = 0)
    case 'y_axis'
        % want a positive y axis crossing
        if abs(t) > 1.e-1 && x(2) > 0 % wait a short time before checking to avoid initial crossing
            isterminal = 1;
        else
            isterminal = 0;
        end
        value=x(1);  % stop at y axis cross (x is 0 or x(1) = 0)
    case 'z_axis'
        % want a positive x axis crossing
        if abs(t) > 1.e-1 && abs(x(1)) > 0 % wait a short time before checking to avoid initial crossing
            isterminal = 1;
        else
            isterminal = 0;
        end
        value=x(1);  % stop at y axis cross (x is 0 or x(1) = 0)
end

end