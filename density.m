function [rho] = density(x)
    delta_x = 10.;
    channel_width = 7.5;
    rho = zeros(size(x));

    [x_sorted, I] = sort(x);
    for i=1:numel(x)
        count = 0;
        for j=(i+1):numel(x)
            if x_sorted(j) < x_sorted(i) + delta_x
                count = count + 1;
            else
                break
            end
        end
        rho(I(i)) = count/(delta_x * channel_width);
    end
end
