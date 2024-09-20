force_a_angle = deg2rad(ShaftProperties.FORCE_A_ANGLE);
force_a = ShaftProperties.FORCE_A_MAGNITUDE * [ ...
    0.0;                    % x
    -sin(force_a_angle);    % y
    -cos(force_a_angle);    % z
];

radius_a = ShaftProperties.DIAMETER_A * 0.5;
torque_a = cross([0.0; radius_a; 0.0], force_a);

radius_b = ShaftProperties.DIAMETER_B * 0.5;






