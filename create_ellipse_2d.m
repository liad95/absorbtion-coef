function [ellipse] = create_ellipse_2d(rows, cols, center_x, center_y, radius_x, radius_y, rotation_angle)


% Define the size of the matrix

% Create a meshgrid for x and y coordinates
[x, y] = meshgrid(1:cols, 1:rows);



% Rotate the coordinates
x_rotated = (x - center_x) * cosd(rotation_angle) - (y - center_y) * sind(rotation_angle) + center_x;
y_rotated = (x - center_x) * sind(rotation_angle) + (y - center_y) * cosd(rotation_angle) + center_y;

% Calculate the equation of the ellipse
ellipse = ((x_rotated - center_x) / radius_x).^2 + ((y_rotated - center_y) / radius_y).^2 <= 1;

% Display the ellipse
figure;
imshow(ellipse);
end
