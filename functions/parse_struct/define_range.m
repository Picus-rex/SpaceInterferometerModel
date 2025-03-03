function [theta_range, theta_x, theta_y] = define_range(angular_extension)
%DEFINE_RANGE

if ~isstruct(angular_extension)
    theta_range = linspace(angular_extension{1}, angular_extension{2}, ...
        angular_extension{3});    
    [theta_x, theta_y] = meshgrid(theta_range, theta_range);
else
    xs = angular_extension.theta_x;
    ys = angular_extension.theta_y;
    theta_range = linspace(xs{1}, xs{2}, xs{3});
    theta_range_y = linspace(ys{1}, ys{2}, ys{3});

    if ys{3} < 2
        error("A minimum of two points is needed. Please consider using y array of 2 elements of 1 order smaller than x array.")
    end

    [theta_x, theta_y] = meshgrid(theta_range, theta_range_y);
end

end