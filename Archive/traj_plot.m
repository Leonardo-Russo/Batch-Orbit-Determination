function traj_plot(t, sol)
% just an additional function to plot the trajectory

    x_ev = sol(2:end, 1);
    y_ev = sol(2:end, 2);

    epochs = t(2:end);
    
    hold on
    
    Rv = 6052;      % km
    r_ev = sqrt(x_ev.^2+y_ev.^2)-Rv;
    
    %plot(epochs, x_ev, 'color', 'green')
    %plot(epochs, y_ev, 'color', 'blue')
    %plot(epochs, r_ev, 'color', 'red')
    plot(x_ev, y_ev, 'color', 'green');
    axis([-9000, 9000, -9000, 9000]);

end