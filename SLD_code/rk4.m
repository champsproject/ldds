

function [ B ] = rk4( t, X, Y, htemp, Delta_brow )

    pos_brow = round(2*t/htemp + 1);
    
    k_1 = h_function(t,           X, Y, htemp, Delta_brow(pos_brow,:));
    k_2 = h_function(t + htemp/2, X + k_1(:,1)/2, Y + k_1(:,2)/2, htemp, Delta_brow(pos_brow+1,:));
    k_3 = h_function(t + htemp/2, X + k_2(:,1)/2, Y + k_2(:,2)/2, htemp, Delta_brow(pos_brow+1,:));
    k_4 = h_function(t + htemp,   X + k_3(:,1), Y + k_3(:,2), htemp, Delta_brow(pos_brow+2,:));
    
    A(:,1) = X;
    A(:,2) = Y;
    
    B = A + (1/6) * (k_1 + 2 * ( k_2 + k_3 ) + k_4);

end

