function [hat_X,hat_P] = KF(x0,P0,sys,Y)
F = sys.F;
G = sys.G;
H = sys.H;
Q = sys.Q;
R = sys.R;
[~,T] = size(Y);
P = P0;
X = x0;
n = length(x0);

hat_X = zeros(n,T);
hat_P = zeros(n,n,T);

for i=1:T
    X = F*X;
    Z_ = H*X;

    P = F*P*F' + G*Q*G';

    K = P*H'*(H*P*H' + R)^-1;
    X = X + K*(Y(:,i) - Z_);
    
    hat_X(:,i) = X;

    P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R*K';
    % P = (eye(n,n) - K*H)*P;
    hat_P(:,:,i) = P;
end
