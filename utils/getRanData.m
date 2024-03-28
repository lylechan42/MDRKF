function [x,y] = getRanData(x0,T,sys,coeff,is_TV)

    x_prev = x0;
    n = size(sys.H,2);
    m = size(sys.H,1);
    x = zeros(n,T);
    y = zeros(m,T);
    
    Delta = 2 * rand - 1;
    F_purt = sys.F + [0, coeff * Delta; 0, 0];
    for t = 1 : T
        x(:,t) = F_purt * x_prev + sys.G * chol(sys.Q) * randn;
        y(:,t) = sys.H * x(:,t) + sys.D * chol(sys.R) * randn(m,1);
        x_prev = x(:,t);
        if is_TV
            Delta = 2 * rand - 1;
            F_purt = sys.F + [0, coeff * Delta; 0, 0];
        end
    end


end

