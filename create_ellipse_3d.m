function [ellipse] = create_ellipse_3d(rows, cols, depth, center_x, center_y, center_z, radius_x, radius_y,radius_z, rotation_angle_x, rotation_angle_y)


% Create a meshgrid for x, y, and z coordinates
[x, y, z] = meshgrid(1:cols, 1:rows, 1:depth);



% Rotate the coordinates around x-axis
x_rotated_x = x;
y_rotated_x = (y - center_y) * cos(rotation_angle_x) - (z - center_z) * sin(rotation_angle_x) + center_y;
z_rotated_x = (y - center_y) * sin(rotation_angle_x) + (z - center_z) * cos(rotation_angle_x) + center_z;

% Rotate the coordinates around y-axis
x_rotated_y = (x_rotated_x - center_x) * cos(rotation_angle_y) + (z_rotated_x - center_z) * sin(rotation_angle_y) + center_x;
y_rotated_y = y_rotated_x;
z_rotated_y = -(x_rotated_x - center_x) * sin(rotation_angle_y) + (z_rotated_x - center_z) * cos(rotation_angle_y) + center_z;

% Calculate the equation of the ellipsoid
ellipsoid_eq = ((x_rotated_y - center_x) / radius_x).^2 + ((y_rotated_y - center_y) / radius_y).^2 + ((z_rotated_y - center_z) / radius_z).^2 <= 1;
ellipsoid = double(ellipsoid_eq);

% Display the ellipsoid
slice(x, y, z, ellipsoid, center_x, center_y, center_z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
colormap(gray);

end