function [rho] = density(x)
    delta_x = 10.;
    rho = zeros(size(x));
    for i=1:numel(x)
        rho(i) = sum((x >= x(i)) .* (x < x(i) + delta_x));
    end
end
