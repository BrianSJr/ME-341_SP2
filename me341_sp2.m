%% Naming Convension:
% Variable names that begin with lower case letter "v" are vectors.
% Variable names that begin with lower case letter "m" are matrices.
% Variable names that begin with neither "v" or "m" are scalars.
% 

% Lengths (L):
Loa = 0.40; % [m] Length from O to A.
Lab = 0.35; % [m] Length from A to B.
Lbc = 0.30; % [m] Length from B to C.
Lob = Loa + Lab; % [m] Length from O to B.
Ls   = Lob + Lbc; % [m] Length from O to C (the whole shaft).

% Diameters (D):
Da  = 0.6; % [m] Diameter of gear A.
Db  = 0.3; % [m] Diameter of gear B.
Ds  = 0.05; % [m] Diameter of shaft.

% Forces (F):
Fa  = 11000.0; % [N] Magnitude of Fa.

% Angles (A)
Aa  = deg2rad(20.0); % [rads] (20 degress) CW from +z-axis about the x-axis.
Ab  = deg2rad(25.0); % [rads] (25 degress) CCW from -z-axis about the x-axis. 

%% Finding Fb:
% To find Fb, we use the sum of the torques about the 
% x-axis:
%
% sum[(T)x] = Da/2 * (Fa)z + Db/2 * (Fb)z = 0
%
% Observing the force on gear A we see:
%               +y
%       Fa       ^ 
%       ▀▄       |
%         ▀▄     |
%           ▀▄ ▄ |
%       Aa ( ▄██ | 
% +z <-----------+
% Therefore:
% >> vFa = -Fa*cos(Aa)j-Fa*sin(Aa)k
% >> (Fa)z = -Fa*sin(Aa)
%
% Observing the force on gear B we see:
%      +y
%       ^      Fb
%       |       ▄▀ 
%       |     ▄▀
%       | ▄ ▄▀
%       | ██▄ ) Ab
% +z <--+-----------> -z
% Therefore:
% >> vFb = -Fb*sin(Ab)j+Fb*cos(Ab)k
% >> (Fb)z = Fb*cos(Ab)
%
% Substituting our (Fa)z and (Fb)z values into the torque equation
% we get:
% sum[(T)x] = Da/2 * (Fa)z + Db/2 * (Fb)z = 0
% <= (Fa)z = -Fa*cos(Aa)
% <= (Fb)z =  Fb*cos(Ab)
% => sum[(T)x] = Da/2 * -Fa * cos(Aa) + Db/2 * Fb * cos(Ab) = 0
%
% Simplify:
% => Da/2 * -Fa * cos(Aa) + Db/2 * Fb * cos(Ab) = 0
% => Db/2 * Fb * cos(Ab) = Da/2 * Fa * cos(Aa)
% => Db * Fb * cos(Ab) = Da * Fa * cos(Aa)
% >> Fb = (Da * Fa * cos(Aa)) / (Db * cos(Ab))
% 
Fb = (Da * Fa * cos(Aa)) / (Db * cos(Ab));

%% Finding Reaction Forces:
% The shaft is in static equilibrium, thus forces Fa and Fb create reaction 
% forces at bearings C and D known as Fc and Fd. To find the Reaction forces 
% draw the FBD in the xy-plane and the yz-plane.

%% FBD xy-plane:
% The FBD in the xy-plane. Note that the arrows are pointing in the direction
% of the sign in the y-direction. They do NOT dictate the sign of the
% component.
% (+y)
%  |
%  |       (Fa)y     (Fb)y
%  |         |         |
%  |         v         v       
%--O---------A---------B---------C--(+x)
%  ^         |         |         ^
%  |<--Loa-->|<--Lab-->|<--Lbc-->|
%(Fo)y                        (Fc)y
%  |
% (-y)

%% FBD yz-plane:
% The FBD in the yz-plane. Note that the arrows are pointing in the direction
% of the sign in the y-direction. They do NOT dictate the sign of the
% component.
% (-z)
%  |
%(Fo)z              (Fb)z
%  |                   |
%  v                   v       
%--O---------A---------B---------C--(+x)
%  |         ^         |         ^
%  |<--Loa-->|<--Lab-->|<--Lbc-->|
%  |       (Fa)z               (Fc)z
%  |
% (+z)

%% Reaction Forces:
% First we create vFa and vFb. Recall we found vFa and vFb earlier.
% >> vFa = -Fa*cos(Aa)j-Fa*sin(Aa)k
% >> vFb = -Fb*sin(Ab)j+Fb*cos(Ab)k
%
vFa = Fa * [0.0; -sin(Aa); -cos(Aa)];
vFb = Fb * [0.0; -sin(Ab);  cos(Ab)];

% To find the vFc, take the moment about the y and z axis at O.
% 
% Moment about y-axis at O:
% sum[(Mo)y] = Loa * (Fa)z + (Loa + Lab) * (Fb)z + (Loa + Lab + Lac) * (Fc)z = 0
%
% Solve for (Fc)z:
% => (Loa + Lab + Lac) * (Fc)z = -(Loa * (Fa)z + (Loa + Lab) * (Fb)z)
% => (Fc)z = -(Loa * (Fa)z + (Loa + Lab) * (Fb)z) / (Loa + Lab + Lac)
% <= Lob = Loa + Lab
% <= Ls = Loa + Lab + Lbc
% >> (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / L
% 
% Moment about z-axis at O:
% sum[(Mo)z] = Loa * (Fa)y + (Loa + Lab) * (Fb)y + (Loa + Lab + Lac) * (Fc)y = 0
% 
% Solve for (Fc)y:
% => (Loa + Lab + Lac) * (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y)
% => (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y) / (Loa + Lab + Lac)
% <= Lob = Loa + Lab
% <= Ls = Loa + Lab + Lbc
% >> (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / L
%
% Note both can be done in parallel using vectors ((Fc)x can be diregarded because 
% both (Fa)x and (Fb)x are zero):
% <= (Fc)x = 0
% <= (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / L
% <= (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / L
% >> vFc = -(Loa * vFa + Lob * vFb) / L
%
vFc = -(Loa * vFa + Lob * vFb) / Ls;

% To find vFo use the sum equilibrium forces as vectors:
%
% sum[vF] = vFo + vFa + vFb + vFc = 0
% 
% Solve for vFo:
% => vFo + vFa + vFb + vFc = 0
% >> vFo = -(vFa + vFb + vFc)
%
vFo = -(vFa + vFb + vFc);


x = linspace(0, Ls, 100);
mV = V(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);
mM = M(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);

figure(1);
plot(x, mV);
figure(2);
plot(x, mM);

%% Shear Force Function
% Shear force on the shaft can be expressed as the following function:
%
% V(x) = u(x)*vFo + u(x-Loa)*vFa + u(x-Lob)*vFb + u(x-Lob)*vFc
% where u(t) = { t < 0: 1, 0 } (step function)
%
function result = V(position, array_of_forces, array_of_scalar_distances) 
    arguments
        position (1, :) {isnumeric},
        array_of_forces (:, :) {isnumeric},
        array_of_scalar_distances (1, :) {isnumeric}
    end
    result = 0;
    for i = 1:size(array_of_forces, 2)
        result = result + array_of_forces(:, i) * (position > array_of_scalar_distances(i));
    end
end

%% Shear Force Function
% Integrating the shear force function reveals the following:
%  
% M(x) = u(x)*x*vFo + u(x-Loa)*(x-Loa)*vFa + u(x-Lob)*(x-Lob)*vFb +
% u(x-Ls)*(x-Ls)*vFc
% where u(t) = { t < 0: 1, 0 } (step function)
%
function result = M(position, array_of_forces, array_of_scalar_distances) 
    arguments
        position (1, :) {isnumeric},
        array_of_forces (:, :) {isnumeric},
        array_of_scalar_distances (1, :) {isnumeric}
    end
    result = 0;
    for i = 1:size(array_of_forces, 2)
        result = result + array_of_forces(:, i) * ((position - array_of_scalar_distances(i)) .* (position > array_of_scalar_distances(i)));
    end
end












