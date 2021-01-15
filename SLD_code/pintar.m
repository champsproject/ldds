

function pintar( X, Y, Lagrangiano, x, y, a, b, tau, potencia, nfig )

    figure(nfig);
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    %set(gcf,'Position',get(0,'ScreenSize'));

    pcolor( X, Y, Lagrangiano );
    %colormap(jet);
    %hold on;
    %plot(a,b,'m*','LineWidth',1);
    % Plot forward hyperbolic orbit 
    %plot(a(:,1),b(:,1),'m*','LineWidth',1);
    %hold on;
    % Plot backward hyperbolic orbit 
    %plot(a(:,2),b(:,2),'g*','LineWidth',1);
    axis tight;
    axis([x(1) x(end) y(1) y(end)]);    
    %title( ['\tau = ', num2str(tau), ', p = ', num2str(potencia) ] );
    %title( ['Magenta = forw iter; Green = backw iter, \tau = ', num2str(tau), ', p = ', num2str(potencia) ] );
    %title( ['Htemp = ', num2str(htemp), ', Xgrid = ', num2str(xgrid), ', Ygrid = ', num2str(ygrid), ', \tau = ', num2str(tau), ', potencia = ', num2str(potencia) ] );
    
    xlabel('x','FontSize',16);
    ylabel('y','FontSize',16);
    shading flat;

end

