function R = prf_contour(Srf, v, X, Y)
%
% R = prf_contour(Srf, v, X, Y)
%
% Returns the pRF profile for vertex v in Srf as a contour matrix. The
% matrices X and Y contain the X and Y coordinates of the visual space.
%
% If output argument R is defined, the function just returns the matrix. 
% If not, it instead displays the pRF profile as a contour plot.
%
% This function requires a Srf analysed with reverse correlation.
%
% 24/02/2020 - Ensured axes are square now (DSS)
%

Rmap = Srf.Rmaps(:,v);
% Reshape vector into a map
R = reshape(Rmap, sqrt(size(Rmap,1)), sqrt(size(Rmap,1)));

if nargout == 0
    contourf(X, Y, R, 50, 'linestyle', 'none');
    axis square
    set(gca, 'fontsize', 12);
    xlabel('Horizontal coordinate (deg)');
    ylabel('Vertical coordinate (deg)');
    colorbar
    title(['X = ' num2str(Srf.Data(2,v)) ', Y = '  num2str(Srf.Data(3,v)) ', R^2 = ' num2str(Srf.Data(1,v))]);
end
