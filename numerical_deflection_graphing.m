%% Naming Convension:
% Variable names that begin with lower case letter "v" are vectors.
% Variable names that begin with lower case letter "m" are matrices.
% Variable names that begin with no lower case letter or the letter "s" are scalars.
%  

% Youngs Modulus:
E = 200.0E9;

% Lengths (L):
Loa = 0.40; % [m] Length from O to A.
Lab = 0.35; % [m] Length from A to B.
Lbc = 0.30; % [m] Length from B to C.
Lob = Loa + Lab; % [m] Length from O to B.
Ls  = Lob + Lbc; % [m] Length from O to C (the whole shaft).

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
% >> (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / Ls
% 
% Moment about z-axis at O:
% sum[(Mo)z] = Loa * (Fa)y + (Loa + Lab) * (Fb)y + (Loa + Lab + Lac) * (Fc)y = 0
% 
% Solve for (Fc)y:
% => (Loa + Lab + Lac) * (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y)
% => (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y) / (Loa + Lab + Lac)
% <= Lob = Loa + Lab
% <= Ls = Loa + Lab + Lbc
% >> (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / Ls
%
% Note both can be done in parallel using vectors ((Fc)x can be diregarded because 
% both (Fa)x and (Fb)x are zero):
% <= (Fc)x = 0
% <= (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / Ls
% <= (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / Ls
% >> vFc = -(Loa * vFa + Lob * vFb) / Ls
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEGIN STUFF DANIEL CHANGED
%
%   this code can be run on its own and does not need to be combined with
%   other file
%
%   Purpose of this file is to graph the deflections and slope of the
%   graph through numerical integration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0, Ls, 100000); %change the third number here to change the number of numerical steps
mV = V(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);
mM = M(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);
mT = T(x, [vFa, vFb], [Loa, Lob], [Da, Db]);

mM(3, :) = -1*mM(3, :); %flip the mM graph going forward

I = (0.25 * pi) * (0.5 * Ds) ^ 4;
disp(1.0*length(x));
dx = (Ls)/(1.0*length(x));


theta_xy = zeros(length(x), 1);
%theta_0=-P*a*b(L+b)/(6EIL) 
theta_0_xy = (vFa(2)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) + vFb(2)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
disp(theta_0_xy)
theta_xy(1) = theta_0_xy;
theta_xy_running = theta_0_xy; %"running" holds the integral sum of the code

y_running = 0;
y_xy = zeros(length(x),1); %saves the integral sum for each x value

for i = 2:length(x)
    theta_xy_running = theta_xy_running + (mM(2, i)*dx)/(E*I); % integral(M(x)/EI, dx)
    theta_xy(i) = theta_xy_running;

    y_running = y_running + theta_xy_running*dx; % integral(theta(x), dx)
    y_xy(i) = y_running;
end

%same thing but in for zx
theta_zx = zeros(length(x), 1);

theta_0_zx = (-vFa(3)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) - vFb(3)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
disp(theta_0_zx)
theta_zx(1) = theta_0_zx;
theta_zx_running = theta_0_zx;

z_running = 0;
z_zx = zeros(length(x),1);

for i = 2:length(x)
    theta_zx_running = theta_zx_running + (mM(3, i)*dx)/(E*I);
    theta_zx(i) = theta_zx_running;

    z_running = z_running + theta_zx_running*dx;
    z_zx(i) = z_running;
end

%figures in the 100's are graph the numerical step approximation
figure(100);
hold on
title("theta xy");
plot(x, theta_xy);
xaxis = yline(0);
xlabel("x (m)");
ylabel("theta (rad)");
hold off

figure(101);
hold on
title("y deflection");
plot(x, y_xy);
xaxis = yline(0);
xlabel("x (m)");
ylabel("y (m)");
hold off

figure(102);
hold on
title("theta zx");
plot(x, -theta_zx);
xaxis = yline(0);
xlabel("x (m)");
ylabel("theta (rad)");
hold off

figure(103);
hold on
title("z deflection");
plot(x, -z_zx);
xaxis = yline(0);
xlabel("x (m)");
ylabel("z (m)");
hold off


%SOLVE FOR THE MINIMUM SHAFT DIAMETER FOR 0.06deg SLOPE

%this code specifically analyzes POINT C as this is where the slope is
%worst

%solve for the minimum diameter
I = (0.25 * pi) * (0.5 * Ds) ^ 4;
disp(1.0*length(x));
dx = (Ls)/(1.0*length(x));


theta_xy = zeros(length(x), 1);
%theta_0=-P*a*b(L+b)/(6EIL) 
theta_0_xy = (vFa(2)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) + vFb(2)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
disp(theta_0_xy)
theta_xy(1) = theta_0_xy;
theta_xy_running = theta_0_xy; %"running" holds the integral sum of the code

y_running = 0;
y_xy = zeros(length(x),1); %saves the integral sum for each x value

for i = 2:length(x)
    theta_xy_running = theta_xy_running + (mM(2, i)*dx)/(E*I); % integral(M(x)/EI, dx)
    theta_xy(i) = theta_xy_running;

    y_running = y_running + theta_xy_running*dx; % integral(theta(x), dx)
    y_xy(i) = y_running;
end

d_h = 1;
d_l = Ds;

while(abs(d_h-d_l) > 0.00000001)

    m = (d_h + d_l) * 0.5;
    if(m == 0)
        b = 0;
        break;
    end

    %calculate for high bound
    I = (0.25 * pi) * (0.5 * d_l) ^ 4;

    theta_xy_low = zeros(length(x), 1);
    theta_0_xy = (vFa(2)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) + vFb(2)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
    theta_xy_low(1) = theta_0_xy;
    theta_xy_running = theta_0_xy; %"running" holds the integral sum of the code
    y_running = 0;
    y_xy_low = zeros(length(x),1); %saves the integral sum for each x value

    for i = 2:length(x)
        theta_xy_running = theta_xy_running + (mM(2, i)*dx)/(E*I); % integral(M(x)/EI, dx)
        theta_xy_low(i) = theta_xy_running;
    
        y_running = y_running + theta_xy_running*dx; % integral(theta(x), dx)
        y_xy(i) = y_running;
    end

    theta_zx_low = zeros(length(x), 1);
    
    theta_0_zx = (-vFa(3)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) - vFb(3)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
    disp(theta_0_zx)
    theta_zx_low(1) = theta_0_zx;
    theta_zx_running = theta_0_zx;
    
    z_running = 0;
    z_zx = zeros(length(x),1);
    
    for i = 2:length(x)
        theta_zx_running = theta_zx_running + (mM(3, i)*dx)/(E*I);
        theta_zx_low(i) = theta_zx_running;
    
        z_running = z_running + theta_zx_running*dx;
        z_zx(i) = z_running;
    end

    %calculate low bound
    I = (0.25 * pi) * (0.5 * m) ^ 4;

    theta_xy_high = zeros(length(x), 1);
    theta_0_xy = (vFa(2)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) + vFb(2)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
    theta_xy_high(1) = theta_0_xy;
    theta_xy_running = theta_0_xy; %"running" holds the integral sum of the code
    y_xy_high = zeros(length(x),1);

    for i = 2:length(x)
        theta_xy_running = theta_xy_running + (mM(2, i)*dx)/(E*I); % integral(M(x)/EI, dx)
        theta_xy_high(i) = theta_xy_running;
    
        y_running = y_running + theta_xy_running*dx; % integral(theta(x), dx)
        y_xy(i) = y_running;
    end

    theta_zx_high = zeros(length(x), 1);
    
    theta_0_zx = (-vFa(3)*Loa*(Lab+Lbc)*(Ls+Lab+Lbc) - vFb(3)*(Loa+Lab)*Lbc*(Ls+Lbc)) / (6 * E * I * Ls);
    disp(theta_0_zx)
    theta_zx_high(1) = theta_0_zx;
    theta_zx_running = theta_0_zx;
    
    z_running = 0;
    z_zx = zeros(length(x),1);
    
    for i = 2:length(x)
        theta_zx_running = theta_zx_running + (mM(3, i)*dx)/(E*I);
        theta_zx_high(i) = theta_zx_running;
    
        z_running = z_running + theta_zx_running*dx;
        z_zx(i) = z_running;
    end

    % 0.06 deg = 0.001047198 rad
    %use pythagorean theorem to get deflection in both directions 
    if(((theta_xy_low(length(x))^2 + (theta_zx_low(length(x))^2))^0.5 - 0.001047198) * (((theta_xy_high(length(x)))^2 + (theta_zx_high(length(x))^2))^0.5 - 0.001047198) < 0)
        d_h = m;
    else
        d_l = m;
    end
end

disp("Minumum diameter for deflection " + num2str(d_h));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END STUFF DANIEL CHANGED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
hold on
title("Shear Force");
plot( ...
    x, mV(1, :), 'r', ...
    x, mV(2, :), 'g', ...
    x, mV(3, :), 'b' ...
    );
xlabel("x (m)");
ylabel("(V)xyz (N)");
hold off

figure(2);
hold on
title("Bending Moment");
plot( ...
    x, mM(1, :), 'r', ...
    x, mM(2, :), 'g', ...
    x, mM(3, :), 'b' ...
    );
xlabel("x (m)");
ylabel("(M)xyz (Nm)");
hold off

%{
figure(3);
plot( ...
    x, mT(1, :), 'r', ...
    x, mT(2, :), 'g', ...
    x, mT(3, :), 'b' ...
    );

I = (0.25 * pi) * (0.5 * Ds) ^ 4;

md2 = Deflect(2, x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls], E, I);
md3 = Deflect(3, x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls], E, I);


md = md2 + md3;

dt = sqrt(sum(md.^2, 1));

md4 = Deflect_s(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls], E, I);
dt2 = sqrt(sum(md4.^2, 1));

figure(7);
plot(x, dt);

figure(8)
plot(x, dt2);

figure(4);
plot( ...
    x, md2(1, :), 'r', ...
    x, md2(2, :), 'g', ...
    x, md2(3, :), 'b' ...
    );

figure(5);
plot( ...
    x, md3(1, :), 'r', ...
    x, md3(2, :), 'g', ...
    x, md3(3, :), 'b' ...
    );

figure(6);
plot( ...
    x, md(1, :), 'r', ...
    x, md(2, :), 'g', ...
    x, md(3, :), 'b' ...
    );

mw = w(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);

figure(12)
plot( ...
    x, mw(1, :), 'r', ...
    x, mw(2, :), 'g', ...
    x, mw(3, :), 'b' ...
    );
%}


%% Load Force Function
function r = w(pvX, mF, vL)
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric}
    end
    r = 0;   
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * (pvX == vL(i));
    end
end

%% Shear Force Function
% Shear force on the shaft can be expressed as the following function:
%
% V(x) = u(x)*vFo + u(x-Loa)*vFa + u(x-Lob)*vFb + u(x-Lob)*vFc
% where u(t) = { t < 0: 1, 0 } (step function)
% 
% Matlab Implementation:
%   By passing a scalar or row vector of positions and arrays representing 
%   the point forces and their position along the beam, the function can 
%   calculate the Shear Force for each position.
% 
% Arguments:
%   svPos - Scalar or row vector (1 by N matrix) representing the position
%           or positions to calcualte the shear force at.
%   mF    - Array of column vectors which represent forces (looks like a matrix).
%   vL    - Array of lengths representing locations of each point force.
% Requirements:
%   - All arguments must satisfy isnumeric.
%   - The number of forces (column vectors) in mF must be equal to the
%     number of lengths (scalars) in vL.
%
function r = V(pvX, mF, vL) 
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths"); 
    r = 0; 
    for i = 1:size(mF, 2)
         r = r + mF(:, i) * (pvX > vL(i));
    end
end

%% Bending Force Function
% Integrating the shear force function reveals the following:
% 
% M(x) = u(x)*x*vFo + u(x-Loa)*(x-Loa)*vFa + u(x-Lob)*(x-Lob)*vFb +
% u(x-Ls)*(x-Ls)*vFc
% where u(t) = { t < 0: 1, 0 } (step function)
%
% Matlab Implementation:
%   By passing a scalar or row vector of positions and arrays representing 
%   the point forces and their position along the beam, the function can 
%   calculate the Shear Force for each position.
% 
% Arguments:
%   svPos - Scalar or row vector (1 by N matrix) representing the position
%           or positions to calcualte the bending moment at.
%   mF    - Array of column vectors which represent forces (looks like a matrix).
%   vL    - Array of lengths representing locations of each point force.
% Requirements:
%   - All arguments must satisfy isnumeric.
%   - The number of forces (column vectors) in mF must be equal to the
%     number of lengths (scalars) in vL.
%
function r = M(pvX, mF, vL) 
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths"); 
   
    r = 0;
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * ((pvX - vL(i)) .* (pvX > vL(i)));
    end
end

%% Torque Function
% Taking the cross product at each gear radius and force according to
% length yeilds the following:
%  
% T(x) = u(x-Loa) * (Da/2)j x vFa + u(x-Lob) * (Db/2)j x vFb
% where u(t) = { t < 0: 1, 0 } (step function)
%
% Matlab Implementation:
%   By passing a scalar or row vector of positions and arrays representing 
%   the point forces and their position along the beam, the function can 
%   calculate the Shear Force for each position.
% 
% Arguments:
%   svPos - Scalar or row vector (1 by N matrix) representing the position
%           or positions to calcualte the torque at.
%   mF    - Array of column vectors which represent forces (looks like a matrix).
%   vL    - Array of lengths representing locations of each point force.
% Requirements:
%   - All arguments must satisfy isnumeric.
%   - The number of forces (column vectors) in mF must be equal to the
%     number of lengths (scalars) in vL.
%
function r = T(pvX, mF, vL, vD)
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric},
        vD (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of diameters");
    r = 0;
    for i = 1:size(mF, 2)
        r = r + cross([0.0; 0.5 * vD(i); 0.0], mF(:, i)) * (pvX > vL(i));
    end
end

function r = Deflect(index, pvX, mF, vL, Es, Is)
    Ls = vL(end);
    Lof = vL(index);
    mM = M(pvX, mF, vL);
    r = (mM / (6 * Es * Is)) .* (ones(3, 1) * (pvX .^ 2 + Lof^2 - 2 * Ls * (Lof - (pvX > Lof) .* (Lof - pvX))));
end

function r = Deflect_s(pvX, mF, vL, Es, Is) 
    r = 0;
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * (((pvX - vL(i)) .^ 3) .* (pvX > vL(i)));
    end
    r = r * 1 / (6 * Es * Is);
end

%% TODO:
% Angular of Twist function (AoT(L, T...))
% Polar Moment of Inertia (PMoI(...))



% 
%
% Yab= ((F*Lb*x) / (6*E*I*Ls)) * (x^2 + Lb^2 - Ls^2)
% Ybc= ((F*La*(Ls - x)) / (6*E*I*Ls)) * (x^2 + Lb^2 - 2*Ls*x)
% <= M = (F*Lb*x) / Ls) maybe
% 
%







