function [x,y] = getCorrData(x0,T,sys,vars)

    [b1, b2, b3, q1, q11, q12, q13] = deal(vars.b1, vars.b2, vars.b3, vars.q1, vars.q11, vars.q12, vars.q13);
    x_prev = x0;
    n = size(sys.H,2);
    m = size(sys.H,1);
    x = zeros(n,T);
    y = zeros(m,T);
    w_prev = 0;
    for t = 1 : T
        w = chol(q1) * randn;
        x(:,t) = sys.F * x_prev + sys.G * w;
        y(1,t) = sys.H(1,:) * x(:,t) + b1*w_prev + chol(q11) * randn;
        y(2,t) = sys.H(2,:) * x(:,t) + b2*w_prev + chol(q12) * randn;
        y(3,t) = sys.H(3,:) * x(:,t) + b3*w_prev + chol(q13) * randn;
        x_prev = x(:,t);
        w_prev = w;
    end


end

