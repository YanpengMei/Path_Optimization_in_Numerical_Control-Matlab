% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% calculate the sum of the curvature of 2D-segment

function sumCur = A12_sumCur2D(var)
load('A2_parameters.mat','cp','step_xy','knot')
load('A2_TupleSeg.mat','TupleSeg')
load('A2_fittedCP_best.mat','fittedCP_best')


% 1. get the fitted segment and 1st- and 2nd-order derivative
% to make sure the programm works, we will only optimize "y-coordinate" of the control points
cp_var = [fittedCP_best(1,:) ; var];

% FCur             fitted NURBS-curve
FCur = nrbmak(cp_var,knot);

% get 1st and 2nd-order derivative
[D, DD] = nrbderiv(FCur);


% 2. get the the number of sampling points
max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2));


% 3. calculate the curvature of fitted segment
% Du                 1st-order derivative of the sampling points in u-direction
% Duu                2nd-order derivative of the sampling points in u-direction
% ut_Fcur            discrete u in segments for the calculation of the curvature

cur = zeros(1,numP_P);

utFcur = linspace(0,1,numP_P);
[~, temp1, temp2] = nrbdeval(FCur, D, DD, utFcur);
Du = temp1';
Duu = temp2';

for j = 1:numP_P
    cur(1,j) = norm(cross(Du(j,:),Duu(j,:)))/norm(Du(j,:))^3;
end
    sumCur = sum(cur);

end