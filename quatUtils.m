function U = quatUtils()
%QUATUTILS  Quaternion helper functions, returned as a struct of handles.
%   Usage:  Q = quatUtils();  R = Q.toR(q);  qe = Q.mul(a,b);  ...
%   Convention: q = [q0 q1 q2 q3], scalar-first, body->geodetic.

    U.mul     = @qmul;       % quaternion product a*b
    U.inv     = @qinv;       % conjugate / inverse (unit quat)
    U.canon   = @canon;      % force scalar part >= 0
    U.toR     = @quat2R;     % quaternion -> rotation matrix (body->world)
    U.Omega   = @quatOmega;  % 4x3 matrix for qdot = 0.5*Omega(q)*w
    U.dcm_z   = @dcm_z;      % body z-axis expressed in world (3rd col of R)
    U.between = @qBetween;   % shortest-arc quaternion rotating v1 -> v2
    U.yaw     = @yawQuat;    % yaw angle -> quaternion about world z
    U.toEulZYX= @quat2eul_zyx;
end

function r = qmul(a,b)
    r = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
         a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
         a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
         a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
end

function r = qinv(q), r = [q(1); -q(2); -q(3); -q(4)]; end

function r = canon(q), if q(1) < 0, r = -q; else, r = q; end, end

function R = quat2R(q)
    q0 = q(1); qv = q(2:4);
    qx = [0 -qv(3) qv(2); qv(3) 0 -qv(1); -qv(2) qv(1) 0];
    R = eye(3) + 2*q0*qx + 2*(qx*qx);
end

function Om = quatOmega(q)
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    Om = [-q1 -q2 -q3;
           q0 -q3  q2;
           q3  q0 -q1;
          -q2  q1  q0];
end

function z = dcm_z(q)
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    z = [2*(q1*q3+q0*q2); 2*(q2*q3-q0*q1); 1-2*(q1^2+q2^2)];
end

function q = qBetween(v1,v2)
    v1=v1/norm(v1); v2=v2/norm(v2); d=dot(v1,v2); ax=cross(v1,v2);
    if d < -0.99999, q=[0;1;0;0]; return; end
    s = sqrt((1+d)*2); q = [s/2; ax/s]; q = q/norm(q);
end

function q = yawQuat(psi), q = [cos(psi/2); 0; 0; sin(psi/2)]; end

function e = quat2eul_zyx(q)
    q=q/norm(q); q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    e = [atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));
         asin(max(-1,min(1, 2*(q0*q2-q3*q1))));
         atan2(2*(q0*q3+q1*q2), 1-2*(q2^2+q3^2))];
end
