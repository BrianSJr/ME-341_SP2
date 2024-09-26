%% Naming Convension:
% Variable names that begin with lower case letter "v" are vectors.
% Variable names that begin with lower case letter "m" are matrices.
% Variable names that begin with no lower case letter or the letter "s" are scalars.
%  

%% TODO:
% 2. Document remaining functions:
%   - Bending Moment
%   - Torque 
%   - Deflection (Singularity)
%   - Deflection (Superposition)
%   - Slope (Singularity)
%   - Slope (Superposition)
% 3. Update Figures (title, x/y labels, legeneds)
%   - Deflection (Singularity)
%   - Deflection (Superposition)
%   - Slope (Singularity)
%   - Slope (Superposition)
%

%% Program Settings:
DIAGRAM_SAMPLES = 100; % number of points generated along the x-axis.

%% Shaft Properties:
% Youngs Modulus:
E = 200.0; % [GPa]

% Lengths (L):
Loa = 400.0; % [mm] Length from O to A.
Lab = 350.0; % [mm] Length from A to B.
Lbc = 300.0; % [mm] Length from B to C.
Lob = Loa + Lab; % [mm] Length from O to B.
Ls  = Lob + Lbc; % [mm] Length from O to C (the whole shaft).

% Diameters (D):
Da  = 600.0; % [mm] Diameter of gear A.
Db  = 300.0; % [mm] Diameter of gear B.
Ds  = 50.0; % [mm] Diameter of shaft.

% Forces (F):
Fa  = 11.0; % [kN] Magnitude of Fa.

% Angles (A)
Aa  = deg2rad(20.0); % [rads] (20 degress) CW from +z-axis about the x-axis.
Ab  = deg2rad(25.0); % [rads] (25 degress) CCW from -z-axis about the x-axis.

% Moment of Inertia (I)
I = (0.25 * pi) * (0.5 * Ds) ^ 4; % [mm^4]

%% Finding Fb:
% To find Fb, we use the sum of the torques about the 
% x-axis:
%
%
% Σ(T)x = 0.5*Da*(Fa)z + 0.5*Db*(Fb)z = 0
%
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
%
%   >> vFa   = -Fa*sin(Aa)j - Fa*cos(Aa)k
%   >> (Fa)z =              - Fa*cos(Aa)
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
%
%   >> vFb   = -Fb*sin(Ab)j + Fb*cos(Ab)k
%   >> (Fb)z =                Fb*cos(Ab)
%
% Substituting our (Fa)z and (Fb)z values into the torque equation
% we get:
% 
%
% Σ(T)x = 0.5*Da*(Fa)z + 0.5*Db*(Fb)z = 0
%
% <= (Fa)z = -Fa*cos(Aa)
% <= (Fb)z =  Fb*cos(Ab)
%
% => Σ(T)x = 0.5*Da*[-Fa*cos(Aa)] + 0.5*Db*[Fb*cos(Ab)] = 0
%
%
% Simplify and Solve for Fb:
%
%
% =>-0.5 * Da * Fa * cos(Aa) + 0.5 * Db * Fb * cos(Ab) = 0
% =>                           0.5 * Db * Fb * cos(Ab) = 0.5 * Da * Fa * cos(Aa)
% =>                                 Db * Fb * cos(Ab) =       Da * Fa * cos(Aa)
% =>                                      Fb           =      (Da * Fa * cos(Aa)) / (Db * cos(Ab))
%
%   >> Fb = (Da * Fa * cos(Aa)) / (Db * cos(Ab))
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
%
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
%

%% FBD yz-plane:
% The FBD in the yz-plane. Note that the arrows are pointing in the direction
% of the sign in the y-direction. They do NOT dictate the sign of the
% component.
%
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
%

%% Reaction Forces:
% First we create vFa and vFb. Recall we found vFa and vFb earlier.
%
%   >> vFa = -Fa*sin(Aa)j - Fa*cos(Aa)k
%   >> vFb = -Fb*sin(Ab)j + Fb*cos(Ab)k
%
vFa = Fa * [0.0; -sin(Aa); -cos(Aa)];
vFb = Fb * [0.0; -sin(Ab);  cos(Ab)];

% To find the vFc, take the moment about the y and z axis at O.
% Moment about y-axis at O:
%
%
% Σ(Mo)y = Loa * (Fa)z + (Loa + Lab) * (Fb)z + (Loa + Lab + Lac) * (Fc)z = 0
%
%
% Solve for (Fc)z:
%
%
% => Loa * (Fa)z + (Loa + Lab) * (Fb)z + (Loa + Lab + Lac) * (Fc)z = 0
% =>                                     (Loa + Lab + Lac) * (Fc)z = -(Loa * (Fa)z + (Loa + Lab) * (Fb)z
% =>                                                         (Fc)z = -(Loa * (Fa)z + (Loa + Lab) * (Fb)z) / (Loa + Lab + Lac)
%
%   >> (Fc)z = -(Loa * (Fa)z + (Loa + Lab) * (Fb)z) / (Loa + Lab + Lac)
%
% <= Lob = Loa + Lab
% <= Ls  = Loa + Lab + Lbc
%
%   >> (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / Ls
% 
% Moment about z-axis at O:
%
%
% Σ(Mo)z = Loa * (Fa)y + (Loa + Lab) * (Fb)y + (Loa + Lab + Lac) * (Fc)y = 0
% 
%
% Solve for (Fc)y:
%
%
% => Loa * (Fa)y + (Loa + Lab) * (Fb)y + (Loa + Lab + Lac) * (Fc)y = 0
% =>                                     (Loa + Lab + Lac) * (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y)
% =>                                                         (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y) / (Loa + Lab + Lac)
%
%   >> (Fc)y = -(Loa * (Fa)y + (Loa + Lab) * (Fb)y) / (Loa + Lab + Lac)
%
% <= Lob = Loa + Lab
% <= Ls  = Loa + Lab + Lbc
%
%   >> (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / Ls
%
% Note both can be done in parallel using vectors ((Fc)x can be diregarded because 
% both (Fa)x and (Fb)x are zero):
%
%
% <= (Fc)x = 0
% <= (Fc)y = -(Loa * (Fa)y + Lob * (Fb)y) / Ls
% <= (Fc)z = -(Loa * (Fa)z + Lob * (Fb)z) / Ls
%
%   >> vFc = -(Loa * vFa + Lob * vFb) / Ls
%
vFc = -(Loa * vFa + Lob * vFb) / Ls;

% To find vFo use the sum equilibrium forces as vectors:
%
%
% Σ[vF] = vFo + vFa + vFb + vFc = 0
% 
%
% Solve for vFo:
%
%
% => vFo + vFa + vFb + vFc = 0
% => vFo                   = -(vFa + vFb + vFc)
%
%   >> vFo = -(vFa + vFb + vFc)
%
vFo = -(vFa + vFb + vFc);

%% Shear, Bending and Torque Graphs:
% Shear, Bending and Torque are calculated at each interval along the beam.
% See each function for more details.

x = linspace(0, Ls, DIAGRAM_SAMPLES); % points along the x axis.

mV = V(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);             % Shear Forces
mM = M(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls]);             % Bending Moments

mT = T(x, [vFa, vFb], [Loa, Lob], [Da, Db]);                    % Torques

mYsp = Ysp(x, [vFa, vFb], [Loa, Lob], Ls, E, I);                % Deflection (Superposition)
mYs  = Ys(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls], E, I);    % Deflection (Singularity)

mYmsp = sqrt(sum(mYsp .^ 2, 1));
mYms = sqrt(sum(mYs .^ 2, 1));

mSsp = Ssp(x, [vFa, vFb], [Loa, Lob], Ls, E, I);
mSs = Ss(x, [vFo, vFa, vFb, vFc], [0, Loa, Lob, Ls], E, I); 

mSmsp = sqrt(sum(mSsp .^ 2, 1));
mSms = sqrt(sum(mSs .^ 2, 1));

fprintf("Θy1 = %f, Θy2 = %f", mSs(2,1), mSs(2, end));
fprintf("Θz1 = %f, Θz2 = %f", mSs(3,1), mSs(3, end));

figure(1);
tiledlayout(3, 4);
nexttile;
plot( ...
    x, mV(1, :), 'r', ...
    x, mV(2, :), 'g', ...
    x, mV(3, :), 'b' ...
    );
title("Shear Force");
xlabel("x (mm)");
ylabel("(V)xyz (kN)");
legend("(V)x", "(V)y", "(V)z");

nexttile;
plot(x, mV(1, :), 'r');
title("Shear Force (x-axis)");
xlabel("x (mm)");
ylabel("(V)x (kN)");
legend("(V)x");

nexttile;
plot(x, mV(2, :), 'g');
title("Shear Force (y-axis)");
xlabel("x (mm)");
ylabel("(V)y (kN)");
legend("(V)y");

nexttile;
plot(x, mV(3, :), 'b');
title("Shear Force (z-axis)");
xlabel("x (mm)");
ylabel("(V)z (kN)");
legend("(V)z");

nexttile
plot( ...
    x, mM(1, :), 'r', ...
    x, mM(2, :), 'g', ...
    x, mM(3, :), 'b' ...
    );
title("Bending Moment");
xlabel("x (mm)");
ylabel("(M)xyz (kNmm)");
legend("(M)x", "(M)y", "(M)z");

nexttile;
plot(x, mM(1, :), 'r');
title("Bending Moment (x-axis)");
xlabel("x (mm)");
ylabel("(M)x (kNmm)");
legend("(M)x");

nexttile;
plot(x, mM(2, :), 'g');
title("Bending Moment (y-axis)");
xlabel("x (mm)");
ylabel("(M)y (kNmm)");
legend("(M)y");

nexttile;
plot(x, mM(3, :), 'b');
title("Bending Moment (z-axis)");
xlabel("x (mm)");
ylabel("(M)z (kNmm)");
legend("(M)z");

nexttile
plot( ...
    x, mT(1, :), 'r', ...
    x, mT(2, :), 'g', ...
    x, mT(3, :), 'b' ...
    );
title("Torque");
xlabel("x (mm)");
ylabel("(T)xyz (kNmm)");
legend("(T)x", "(T)y", "(T)z");

nexttile
plot(x, mT(1, :), 'r');
title("Torque (x-axis)");
xlabel("x (mm)");
ylabel("(T)x (kNmm)");
legend("(T)x");

nexttile
plot(x, mT(2, :), 'g');
title("Torque (y-axis)");
xlabel("x (mm)");
ylabel("(T)y (kNmm)");
legend("(T)y");

nexttile
plot(x, mT(3, :), 'b');
title("Torque (z-axis)");
xlabel("x (mm)");
ylabel("(T)z (kNmm)");
legend("(Y)z");

figure;
tiledlayout(4, 4);
nexttile;
plot( ...
   x, mYsp(1, :), 'r', ...
   x, mYsp(2, :), 'g', ...
   x, mYsp(3, :), 'b' ...
   );
title("Deflection (Super Position)");
xlabel("x (mm)");
ylabel("(Y)xyz (mm)");
legend("(Y)x","(Y)y","(Y)z");

nexttile;
plot(x, mYsp(1, :), 'r');
title("Deflection (Super Position) (x-axis)");
xlabel("x (mm)");
ylabel("(Y)x (mm)");
legend("(Y)x");

nexttile;
plot(x, mYsp(2, :), 'g');
title("Deflection (Super Position) (y-axis)");
xlabel("x (mm)");
ylabel("(Y)y (mm)");
legend("(Y)y");

nexttile;
plot(x, mYsp(3, :), 'b');
title("Deflection (Super Position) (z-axis)");
xlabel("x (mm)");
ylabel("(Y)z (mm)");
legend("(Y)z");

nexttile;
plot( ...
   x, mYs(1, :), 'r', ...
   x, mYs(2, :), 'g', ...
   x, mYs(3, :), 'b' ...
   );
title("Deflection (Singularity)");
xlabel("x (mm)");
ylabel("(Y)xyz (mm)");
legend("(Y)x","(Y)y","(Y)z");

nexttile;
plot(x, mYs(1, :), 'r');
title("Deflection (Singularity) (x-axis)");
xlabel("x (mm)");
ylabel("(Y)x (mm)");
legend("(Y)x");

nexttile;
plot(x, mYs(2, :), 'g');
title("Deflection (Singularity) (y-axis)");
xlabel("x (mm)");
ylabel("(Y)y (mm)");
legend("(Y)y");

nexttile;
plot(x, mYs(3, :), 'b');
title("Deflection (Singularity) (z-axis)");
xlabel("x (mm)");
ylabel("(Y)z (mm)");
legend("(Y)z");

nexttile;
plot( ...
   x, mSs(1, :), 'r', ...
   x, mSs(2, :), 'g', ...
   x, mSs(3, :), 'b' ...
   );
title("Slope (Singularity)");
xlabel("x (mm)");
ylabel("(S)xyz (mm)");
legend("(S)x","(S)y","(S)z");

nexttile;
plot(x, mSs(1, :), 'r');
title("Slope (Singularity) (x-axis)");
xlabel("x (mm)");
ylabel("(S)x (mm)");
legend("(S)x");

nexttile;
plot(x, mSs(2, :), 'g');
title("Slope (Singularity) (y-axis)");
xlabel("x (mm)");
ylabel("(S)y (mm)");
legend("(S)y");

nexttile;
plot(x, mSs(3, :), 'b');
title("Slope (Singularity) (z-axis)");
xlabel("x (mm)");
ylabel("(S)z (mm)");
legend("(S)z");

nexttile;
plot( ...
   x, mSsp(1, :), 'r', ...
   x, mSsp(2, :), 'g', ...
   x, mSsp(3, :), 'b' ...
   );
title("Slope (Singularity)");
xlabel("x (mm)");
ylabel("(S)xyz (rads)");
legend("(S)x","(S)y","(S)z");

nexttile;
plot(x, mSsp(1, :), 'r');
title("Slope (Singularity) (x-axis)");
xlabel("x (mm)");
ylabel("(S)x (rads)");
legend("(S)z");

nexttile;
plot(x, mSsp(2, :), 'g');
title("Slope (Singularity) (y-axis)");
xlabel("x (mm)");
ylabel("(S)y (rads)");
legend("(S)y");

nexttile;
plot(x, mSsp(3, :), 'b');
title("Slope (Singularity) (z-axis)");
xlabel("x (mm)");
ylabel("(S)z (rads)");
legend("(S)z");

figure;
tiledlayout(2, 2);
nexttile;
plot(x, mYmsp, 'k');
title("Deflection (Super Position) (Magnitude)");
xlabel("x (mm)");
ylabel("Ysp (mm)");
legend("|Ysp|");

nexttile;
plot(x, mYmsp, 'k');
title("Deflection (Singularity) (Magnitude)");
xlabel("x (mm)");
ylabel("Ys (mm)");
legend("|Ys|");

nexttile;
plot(x, mYmsp, 'k');
title("Slope (Super Position) (Magnitude)");
xlabel("x (mm)");
ylabel("Ssp (rads)");
legend("|Ssp|");

nexttile;
plot(x, mYmsp, 'k');
title("Slope (Singularity) (Magnitude)");
xlabel("x (mm)");
ylabel("Ss (rads)");
legend("|Ss|");



% nexttile;
% plot(x, mSsp(3, :), 'b');
% title("Slope (Singularity) (Magnitude)");
% xlabel("x (mm)");
% ylabel("(S) (mm)");
% legend("|S|");

% vv update/create plots (title, x/y labels, legends) vv
% figure(5);
% plot( ...
%    x, mYs(1, :), 'r', ...
%    x, mYs(2, :), 'g', ...
%    x, mYs(3, :), 'b' ...
%    );
% figure(6);
% plot( ...
%    x, mSsp(1, :), 'r', ...
%    x, mSsp(2, :), 'g', ...
%    x, mSsp(3, :), 'b' ...
%    );
% figure(7);
% plot( ...
%    x, mSs(1, :), 'r', ...
%    x, mSs(2, :), 'g', ...
%    x, mSs(3, :), 'b' ...
%    );

%% Shear Force Function
% Explanation:
% Using singularity functions we can express our load function as:
%
%
% vw(x) = vFo*<x>^(-1) + vFa*<x-Loa>^(-1) + vFb*<x-Lob>^(-1) + vFc*<x-Ls>^(-1).
% 
%
% Singularity Functions:
%
%
% <= <x-a>ⁿ = { δ{|n+1|}(x-a)  when n <  0 } 
%             { (x-a)ⁿH(x-a)   when n >= 0 }
% Where δ{i}(x) is the i-th derivative of the dirac delta function. 
% (eg. i=0 --> δ(x), i=1 --> δ'(x), i=2 --> δ''(x), ...). Note |n+1| in the 
% formula when applying it.     
%
% <= δ(x)   = { 0   when x <  0 }
%             { 1   when x >= 0 }
%           
% <= H(x)   = { 1   when x == 0 }
%             { 0   when x != 0 } 
%
% Integrating vw(x) yeilds the following:
% 
%         ⌠             ⌠                ⌠                    ⌠                    ⌠
% vV(x) = ⌡vw(x)dx = vFo⌡<x>^(-1)dx + vFa⌡<x-Loa>^(-1)dx + vFb⌡<x-Lob>^(-1)dx + vFc⌡<x-Ls>^(-1)dx.
%
%    ⌠           { < x - a >^(n+1)          when n <  0 }
% <= ⌡<x-a>ⁿdx = { < x - a >^(n+1)/(n+1)    when n >= 0 }
%
%
%      ⌠
%   >> ⌡(<x-a>^(-1))dx = <x-a>^0 + C
%
% => vV(x) = vFo<x>^0 + vFa<x-Loa>^0 + vFb<x-Lob>^0 + vFc<x-Ls>^0 + vCv
% => vV(x) = vFo*H(x) + vFa*H(x-Loa) + vFb*H(x-Lob) + vFc*H(x-Ls) + vCv
%
%   >> V(x) = vFo*H(x) + vFa*H(x-Loa) + vFb*H(x-Lob) + vFc*H(x-Ls) + vCv
%
% Solve for when V(0) = 0i + 0j + 0k:
% 
%
% => vV(0) = vFo*H(0) + vFa*H(0-Loa) + vFb*H(0-Lob) + vFc*H(0-Ls) + vCv = 0i + 0j + 0k
%
% <= H(0)     = 0
% <= H(0-Loa) = H(-Loa) = 0
% <= H(0-Lob) = H(-Lob) = 0
% <= H(0-Ls)  = H(0)    = 0
%
% => vV(0) = vFo*(0) + vFa*(0) + vFb*(0) + vFc*(0) + vCv = 0i + 0j + 0k
% 
%   >> vCv   = 0i + 0j + 0k
%   >> vV(x) = vFo*H(x) + vFa*H(x-Loa) + vFb*H(x-Lob) + vFc*H(x-Ls)
%
% Matlab Implementation:
%   By passing a scalar or row vector of positions and arrays representing 
%   the point forces and their position along the beam, the function can 
%   calculate the Shear Force for each position.
% 
% Arguments:
%   svX - Scalar or row vector (1 by N matrix) representing the position
%         or positions to calcualte the shear force at.
%   mF  - Array of column vectors which represent forces (looks like a matrix).
%   vL  - Array of lengths representing locations of each point force.
% Requirements:
%   - All arguments must satisfy isnumeric.
%   - The number of forces (column vectors) in mF must be equal to the
%     number of lengths (scalars) in vL.
%
function r = V(svX, mF, vL) 
    arguments
        svX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths"); 
    r = 0; 
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * (svX > vL(i));
    end
end


%% Bending Force Function:
% Recall our shear force function vV(x) with the singularity function definition:
%
%   >> vV(x) = vFo<x>^0 + vFa<x-Loa>^0 + vFb<x-Lob>^0 + vFc<x-Ls>^0 + vCv
% 
% Integrating the shear force function yeilds the following:
% Note that we use the cross product between the arm vector and the force
% at that position. 
%
%                ⌠                  ⌠                ⌠                    ⌠                    ⌠
% => vM(x) = i x ⌡vV(x)dx = i x (vFo⌡<x>^(-1)dx + vFa⌡<x-Loa>^(-1)dx + vFb⌡<x-Lob>^(-1)dx + vFc⌡<x-Ls>^(-1)dx)
%    
%      ⌠
%   >> ⌡(<x-a>^0)dx = <x-a>^1 + C
%
% => vM(x) = (<x>^1)i x vFo + (<x-Loa>^1)i x vFa + (<x-Lob>^1)i x vFb + (<x-Ls>^1)i x vFc 
% => vM(x) = (xH(x))i x vFo + (x-Loa)H(x-Loa)i x vFa + (x-Lob)H(x-Loa)i + vFb + (x-Ls)H(x-Ls)i x vFc + vCm
% 
%   >> vM(x) = (xH(x))i x vFo + (x-Loa)H(x-Loa)i x vFa + (x-Lob)H(x-Loa)i + vFb + (x-Ls)H(x-Ls)i x vFc + vCm
%  
% Solve for vCm when vM(0) = 0i + 0j + 0k:
%
% => vM(0) = 0i + 0j + 0k + vCm
%
%   >> vCM = 0i + 0j + 0k
%   >> vM(x) = (xH(x))i x vFo + (x-Loa)H(x-Loa)i x vFa + (x-Lob)H(x-Loa)i + vFb + (x-Ls)H(x-Ls)i x vFc
%
% Matlab Implementation:
%   By passing a scalar or row vector of positions and arrays representing 
%   the point forces and their position along the beam, the function can 
%   calculate the Bending Moment Force for each position.
% 
% Arguments:
%   svPos - Scalar or row vector (1 by N matrix) representing the position
%           or positions to calcualte the Bending Moment at.
%   mF    - Array of column vectors which represent forces (looks like a matrix).
%   vL    - Array of lengths representing locations of each point force.
%
% Requirements:
%   - All arguments must satisfy isnumeric.
%   - The number of forces (column vectors) in mF must be equal to the
%     number of lengths (scalars) in vL and the number of rows must be
%     either 1 or 3.
%
function r = M(pvX, mF, vL) 
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths"); 
    assert(size(mF, 1) == 1 || size(mF, 1) == 3, "must be 3D or 1D");
    r = 0;
    if(size(mF, 1))
        for i = 1:size(mF, 2)
            r = r + mF(:,i) * ((pvX - vL(i)) .* (pvX > vL(i)));
        end
    else
        for i = 1:size(mF, 2)
            r = r + transpose(cross([1.0, 0.0, 0.0], mF(:, i))) * ((pvX - vL(i)) .* (pvX > vL(i)));
        end
    end
end

%% Torque Function vv needs doc update vv
% Taking the cross product at each gear radius and force according to
% length yeilds the following:
%  
% T(x) = (Da/2)H(x-Loa)j x vFa + (Db/2)H(x-Lob)j x vFb
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
        mF (3, :) {isnumeric},
        vL (1, :) {isnumeric},
        vD (1, :) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    assert(size(mF, 1) == 1 || size(mF, 1) == 3, "must be 3D or 1D");
    r = 0;
    for i = 1:size(mF, 2)
        r = r + cross([0.0; 0.5 * vD(i); 0.0], mF(:, i)) * (pvX > vL(i));
    end
end

%% Delfection (Method of Superposition):
% vv needs docs vv
function r = Ysp(pvX, mF, vL, Ls, Es, Is)
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric},
        Ls (1, 1) {isnumeric},
        Es (1, 1) {isnumeric},
        Is (1, 1) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    r = 0;
    pvXsqrd = pvX .^ 2;
    Lssqrd = Ls ^ 2;
    for i = 1:size(mF, 2)
        Hx = pvX > vL(i);
        Lof = vL(i);
        Lfc = Ls - Lof; 
        r = r - mF(:, i) * ((~Hx) .* pvX .* (Lfc * (pvXsqrd + Lfc ^ 2 - Lssqrd)) + Hx .* (Ls - pvX) .* (Lof * (pvXsqrd + Lof ^ 2 - 2 * Ls * pvX)));
    end
    r = r * 1 / (6 * Es * Is * Ls);
end

%% Deflection (Method of Singularity):
%
function r = Ys(pvX, mF, vL, Es, Is) 
   arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric},
        Es (1, 1) {isnumeric},
        Is (1, 1) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    r = 0;
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * (((pvX - vL(i)) .^ 3) .* (pvX > vL(i)));
    end
    r = r - (mF(:, 1:(end - 1)) * transpose((vL(end) - vL(1:(end - 1))) .^ 3)) * pvX / vL(end);
    r = r * 1 / (6 * Es * Is);  
end 

%% Slope Function (Method of Superposition)
% vv not implemented (duh) vv
function r = Ssp(pvX, mF, vL, Ls, Es, Is)
    arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric},
        Ls (1, 1) {isnumeric},
        Es (1, 1) {isnumeric},
        Is (1, 1) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    r = 0;
    pvXsqrd = pvX .^ 2;
    Lssqrd = Ls ^ 2;
    for i = 1:size(mF, 2)
        Hx = pvX > vL(i);
        Lof = vL(i);
        Lfc = Ls - Lof; 
        r = r - mF(:, i) * ((~Hx) .* (Lfc * (3 * pvXsqrd + Lfc ^ 2 - Lssqrd)) - Hx .* (Lof * (3 * pvXsqrd + Lof ^ 2 + 2 * Lssqrd - 6 * Ls * pvX)));
    end
    r = r * 1 / (6 * Es * Is * Ls);
end

%% Slope Function (Method of Singularity) 
% vv works vv need docs vv
function r = Ss(pvX, mF, vL, Es, Is) 
   arguments
        pvX (1, :) {isnumeric},
        mF (:, :) {isnumeric},
        vL (1, :) {isnumeric},
        Es (1, 1) {isnumeric},
        Is (1, 1) {isnumeric}
    end
    assert(size(mF, 2) ~= size(vL, 1), "the number of forces must be equal to the number of lengths");
    r = 0;
    for i = 1:size(mF, 2)
        r = r + mF(:, i) * (((pvX - vL(i)) .^ 2) .* (pvX > vL(i)));
    end
    r = r - (mF(:, 1:(end - 1)) * transpose((vL(end) - vL(1:(end - 1))) .^ 3)) / (3 * vL(end));
    r = r * 1 / (2 * Es * Is);  
end 










