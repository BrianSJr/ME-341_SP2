%% class ShaftProperties:
% Stores constant data relavent to the shaft according to the figure.
% Can be used accross multiple files.
classdef ShaftProperties
    properties(Constant)
        LENGTH_OA           = 0.40; % m
        LENGTH_AB           = 0.35; % m
        LENGTH_BC           = 0.30; % m
        LENGTH              = 1.05; % m (taken from summing LENGTH_OA, LENGTH_AB, and LENGTH_BC)
        DIAMETER_A          = 0.6; % m
        DIAMETER_B          = 0.3; % m
        DIAMETER            = 0.05; % m
        FORCE_A_MAGNITUDE   = 11000.0; % N
        FORCE_A_ANGLE       = 20.0; % degs from +z-axis about the x-axis (polarity of force is pointing the opposite direction with respect to the angle)
        FORCE_B_ANGLE       = 25.0; % degs from -z-axis about the x-axis (polarity of force is pointing the opposite direction with respect to the angle)
    end
end