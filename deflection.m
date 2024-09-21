%This is the youngs modulus
E = 200*10^9;

%Moment of inertia for circle
I= (pi/4)*(Shaft_diameter/2)^4


%This is the function for deflection. hasn't been adapted to use vectors
%and matrices

%The formula is taken from the book

%axis_distance is distance from the left of the beam, and point_distance is
%where the force is

%The reason that there is TBD is because I haven't when x goes past where
%the force is. That part of the expression is different

function deflection=calculateDeflection(BendingMoment,E,axis_distance,I,point_distance,TBD)
    deflection = (BendingMoment /(6*E*I))*(axis_distance^2+ point_distance^2+TBD)
end