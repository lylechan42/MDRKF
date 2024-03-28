function [sys,x0,P0] = getSys(dt,choice,vars)

if nargin < 3
    vars = 0; % You can change this default value as needed
end

switch choice
    case 0
        sys.F = 0.95*[1 dt;
                      0 1];
        sys.G = [1 0;
                 0 1];
        sys.H = [1 0;
                 1 0;
                 0 1];
        sys.D = 1;
        sys.Q = 0.01*[dt^3/3 dt^2/2;
                      dt^2/2 dt];
        sys.R = [.25    .125   .00125;
                 .125   .25    .00125;
                 .00125 .00125 .0025];
        x0 = [10;0.5];
        P0 = [1    1/dt;
              1/dt 2/dt^2];
    case 1
        rho = vars;
        sys.F = [1 dt;
                 0 1];
        sys.G = [dt^2/2;
                 dt];
        sys.H = [1 0;
                 1 0];
        sys.D = 1;
        sys.Q = 10;
        sys.R = [10              rho*sqrt(10*9);
                 rho*sqrt(10*9)  9];
        x0 = [10;5];
        P0 = [100 0;
              0   25];
    case 2
        rho = vars;
        sys.F = [1 dt;
                 0 1];
        sys.G = [dt^2/2;
                 dt];
        sys.H = [1 0;
                 0 1];
        sys.D = 1;
        sys.Q = 10;
        sys.R = [10              rho*sqrt(10*2);
                 rho*sqrt(10*2)  2];
        x0 = [10;5];
        P0 = [100 0;
              0   25]; 
    case 3
        % q1:sigma_xi; q11:sigma_eta1;
        % q12:sigma_eta2; q13:sigma_eta3;
        [b1, b2, b3, q1, q11, q12, q13] = deal(vars.b1, vars.b2, vars.b3, vars.q1, vars.q11, vars.q12, vars.q13);

        sys.F = [1 dt dt^2/2;
                 0 1  dt;
                 0 0  1];
        sys.G = [dt^2/2; dt; 1];
        sys.H = [1 0 0;
                 0 1 0;
                 0 0 1];
        sys.D = 1;
        sys.Q = 0.01;
%         sys.R = [b1^2*q1+q11, b1*b2*q1,    b1*b3*q1;
%                  b2*b1*q1,    b2^2*q1+q12, b2*b3*q1;
%                  b3*b1*q1,    b3*b2*q1,    b3^2*q1+q13];
        sys.R = diag([q11,q12,q13]);
        x0 = [0;0;0];
        P0 = eye(3);
end

end

