function [x,y] = getData(x0,T,sys)

    x_prev = x0;
    n = size(sys.H,2);
    m = size(sys.H,1);
    x = zeros(n,T);
    y = zeros(m,T);

    for t = 1 : T
        x(:,t) = sys.F * x_prev + sys.G * chol(sys.Q) * randn(n,1);
        y(:,t) = sys.H * x(:,t) + sys.D * chol(sys.R) * randn(m,1);
        x_prev = x(:,t);
    end


end

