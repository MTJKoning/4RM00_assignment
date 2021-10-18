function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX Cmu Ti
% variables
global y u v T m_in m_out y_v F_u k eps x_u F_v p

% Set velocity at the inlet
for J = 1:NPJ+2
    u(2,1:NPJ+2) = U_IN; % inlet
    %u(2,J) = U_IN*1.5*(1.-(2.*(y(J)-YMAX/2.)/YMAX)^2); % inlet
end

% Set k and eps at the inlet
k(1,1:NPJ+2)     = 1.5*(U_IN*Ti)^2; % at inlet
eps(1,1:NPJ+2)   = Cmu^0.75 *k(1,1:NPJ+2).^1.5/(0.07*YMAX*0.5); % at inlet

% Fix temperature at the walls in Kelvin
T(1:NPI+2,1) = 283.; % bottom wall
T(1:NPI+2,NPJ+2) = 283.; % top wall

% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in = 0.;
m_out_u = 0.;
m_out_v = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out_u = m_out_u + F_u(NPI+1,J)*AREAw;
end

for I = 2:NPI+1
    i=I;
    AREAs = x_u(i+1) - x_u(i);
    m_out_v = m_out_v + F_v(I,NPJ+1)*AREAs;
end
    
% end: globcont()==========================================================

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
% Left side no wall
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1)*m_in/m_out_u;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
% Top side no wall
u(2:NPI+1,NPJ+1) = u(2:NPI+1,NPJ+1);
v(2:NPI+1,NPJ+1) = v(2:NPI+1,NPJ+1); 

% Pressure left side
p(NPI+2,2:NPJ+1) = p(NPI+1,2:NPJ+1);
p(NPI+2,2:NPJ+1) = p(NPI+1,2:NPJ+1);
% Top side no wall
p(2:NPI+1,NPJ+2) = p(2:NPI+1,NPJ+1);
p(2:NPI+1,NPJ+2) = p(2:NPI+1,NPJ+1);
% % Pressure right side
% p(1,2:NPJ+1) = p(2,2:NPJ+1);
% p(1,2:NPJ+1) = p(2,2:NPJ+1);
% % down side no wall
% p(2:NPI+1,1) = p(2:NPI+1,2);
% p(2:NPI+1,1) = p(2:NPI+1,2);


k(NPI+2,2:NPJ+1) = k(NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
T(NPI+2,1:NPJ+2) = T(NPI+1,1:NPJ+2);
end
