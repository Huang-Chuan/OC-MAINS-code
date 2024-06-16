function [] = plot_trajectory_MAINS_OCMAINS(mains_out, ocmains_out, in_data, magdata)

    

    X = mains_out.x_h(:, 1:end)';
    X_ocmains = ocmains_out.x_h(:, 1:end)';

    Ps = mains_out.diag_P(:, 1:end).';
    ref.rotations = in_data.gt.ori(:, :, 1:end);
    ref.positions = in_data.gt.pos(:, 1:end);
    %nis      = mains_out.nis(1:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%               position plot             %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[100 100 800 600]);
    fontsize(18,"points")
    plot(X(:,1),X(:,2), 'k', 'LineWidth', 1.4);
    

    hold on
    plot(X_ocmains(:, 1), X_ocmains(:,2), 'color','#EDB120', 'LineWidth', 1.4);
    if exist('ref','var')
        plot(ref.positions(1,:),ref.positions(2,:), 'm', 'LineWidth', 1.4);
    end



    
    s=pcolor(magdata.X_grid, magdata.Y_grid, reshape(magdata.mag_norm, size(magdata.X_grid)));
    s.EdgeColor= 'none'; 
    s.AlphaData= 0.1./reshape(magdata.ysd, size(magdata.X_grid)); 
    s.FaceAlpha='flat';
    % create a color bar
    c= colorbar('eastoutside', 'FontSize', 14);
    c.Label.String = 'Magnetic field magnitude [\mu T]';
    plot3(X(1,1),X(1,2), 0.8,  'rs', 'LineWidth', 2);
    %plot3(X(end,1),X(end,2), 0.8,'ro','LineWidth', 2);



    ylim([-1.5 1.5])
    grid on;
    title('Horizontal plane trajectory','FontName','Times New Roman','FontSize',18)
    legend('MAINS trajectory', 'OC-MAINS trajectory', 'True trajectory','','Start point',...
        'FontName','Times New Roman', 'location', 'best','FontSize',10)
    xlabel('x [m]','FontName','Times New Roman','FontSize',14)
    ylabel('y [m]','FontName','Times New Roman','FontSize',14)

   
 
end