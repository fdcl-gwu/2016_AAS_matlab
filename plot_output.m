function plot_output(sol_output,reach_struct,min_reach)

% generate and save the handles to a bunch of figures
num_figs = 6;
fontsize = 18;
fontname = 'Times';

fig_title = {'Poincar\`e Section','Poincar\`e Section', 'Reach Distance','Reach Distance','Reach Distance','Trajectory'};
fig_xlabel = {'$x$','$z$','$\phi_1$','$\phi_2$','$\phi_3$','$x$'};
fig_ylabel = {'$\dot x$','$\dot z$','dist','dist','dist','$y$'};

fig_handle = zeros(num_figs);

for ii = 1:num_figs
    fig_handle(ii) = figure();
    set(0, 'CurrentFigure', fig_handle(ii)) 
    hold all
    grid on
    title(fig_title(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    xlabel(fig_xlabel(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    ylabel(fig_ylabel(ii),'interpreter','latex','FontName',fontname,'FontSize',fontsize);
    set(gca,'FontName',fontname,'FontSize',fontsize);
end

reach_state = cat(1,reach_struct(:).reach_end);
pmap = sol_output(1).constants.pmap;

min_state = min_reach.xf;

% Poincare section x vs xdot
set(0,'CurrentFigure',fig_handle(1));
plot(reach_state(:,pmap(1)),reach_state(:,pmap(3)),'rx')

% Poincare section z vs zdot
set(0,'CurrentFigure',fig_handle(2));
plot(reach_state(:,pmap(2)),reach_state(:,pmap(4)),'rx')

% Distance vs Phi1
phi_d = cat(1,reach_struct(:).phi_d)*180/pi;
dist = [reach_struct(:).dist];
min_phi_d = reach_struct(min_reach.index).phi_d*180/pi;
min_dist = min_reach.dist;

set(0,'CurrentFigure',fig_handle(3));
plot(phi_d(:,1),dist,'rx');
plot(min_phi_d(1),min_dist,'go')
% Distance vs Phi2
set(0,'CurrentFigure',fig_handle(4));
plot(phi_d(:,2),dist,'rx');
plot(min_phi_d(2),min_dist,'go')

% Distance vs Phi3
set(0,'CurrentFigure',fig_handle(5));
plot(phi_d(:,3),dist,'rx');
plot(min_phi_d(3),min_dist,'go')

% Trajectory
set(0,'CurrentFigure',fig_handle(6));
vertex_plotter(sol_output(1).constants.F,sol_output(1).constants.V,fig_handle(6));
plot3(min_reach.traj(:,1),min_reach.traj(:,2),min_reach.traj(:,3),'k');

% for ii = 1:length(sol_output)
% %     plot_seg(sol_output(ii).t,sol_output(ii).x_i,sol_output(ii).h_i,sol_output(ii).xm,sol_output(ii).hm, sol_output(ii).x0, sol_output(ii).h0);
%     set(0,'CurrentFigure',poincare1_fig);
%     plot(sol_output(ii).x_i(end,1,end),sol_output(ii).x_i(end,4,end),'ro')
%     
%     set(0,'CurrentFigure',poincare2_fig);
%     plot(sol_output(ii).x_i(end,3,end),sol_output(ii).x_i(end,6,end),'ro')
% end

% plot the constraint satisfaction for each iteration
end

function plot_seg(t,x_i,h_i,xm,hm, x0, h0,constants)
% plot the segments to see if they're correct
traj_fig = figure(100);
grid on
hold all

costate_fig = figure(200);
grid on
hold all
state_fig = figure(300);
grid on
hold all
for ii = 1:size(x_i,3)
    set(0,'CurrentFigure',traj_fig);
    plot3(x_i(:,1,ii),x_i(:,2,ii),x_i(:,3,ii))
    if ii == 1
        plot3(x0(1),x0(2),x0(3),'ko')
        text(x0(1),x0(2),x0(3),sprintf('%s%d%s','x',ii-1,'+'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        statef = [x_i(end,:,ii) h_i(end,:,ii)];
        
    else
        plot3(xm(1,ii-1),xm(2,ii-1),xm(3,ii-1),'ko')
        text(xm(1,ii-1),xm(2,ii-1),xm(3,ii-1),sprintf('%s%d%s','x',ii-1,'+'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        statef = [x_i(end,:,ii) h_i(end,:,ii)];
        
    end
    plot3(statef(1),statef(2),statef(3),'rs','markersize',10)
    text(statef(1),statef(2),statef(3),sprintf('%s%d%s','x',ii,'-'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
    % calculate the x_minus, h_minus terms
    
    % costate
    set(0,'CurrentFigure',costate_fig);
    for ind = 1:size(h_i,2)
        subplot(6,1,ind)
        grid on
        hold all
        plot(t(ii,:),h_i(:,ind,ii))
        title(sprintf('%s %d','lambda',ind));
        
        if ii == 1
            plot(t(ii,1),h0(ind),'ko')
            text(t(ii,1),h0(ind),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        else
            plot(t(ii,1),hm(ind,ii-1),'ko')
            text(t(ii,1),hm(ind,ii-1),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        end
        
        plot(t(ii,end),statef(6+ind),'rs','markersize',10)
        text(t(ii,end),statef(6+ind),sprintf('%s%d%s%d%s','$\lambda_{',ind,'}^{',ii,'}-$'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
    end
    % state
    set(0,'CurrentFigure',state_fig);
    for ind = 1:size(h_i,2)
        subplot(6,1,ind)
        grid on
        hold all
        plot(t(ii,:),x_i(:,ind,ii))
        title(sprintf('%s %d','x',ind));
        
        
        if ii == 1
            plot(t(ii,1),x0(ind),'ko')
            text(t(ii,1),x0(ind),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        else
            plot(t(ii,1),xm(ind,ii-1),'ko')
            text(t(ii,1),xm(ind,ii-1),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii-1,'}+$'),'HorizontalAlignment','left','verticalalignment','top','interpreter','latex')
        end
        plot(t(ii,end),statef(ind),'rs','markersize',10)
        text(t(ii,end),statef(ind),sprintf('%s%d%s%d%s','$x_{',ind,'}^{',ii,'}-$'),'HorizontalAlignment','right','verticalalignment','bottom','interpreter','latex')
    end
    
    
end


end % end of plot segment function

function plot_constraints(sol_output)

end