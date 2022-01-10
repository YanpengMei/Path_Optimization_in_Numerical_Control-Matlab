% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% calculate sum of the curvature of the iso-segment


function sumCur = A10_sumCurIso(cp_var)
load('A2_parameters.mat','knots','step_xy','cp')
load('A2_TupleSeg.mat','TupleSeg')

% 1. get surface
surf = nrbmak(cp_var,knots);


% 2. get 1st and 2nd derivative of the surface
[surfD, surfDD] = nrbderiv(surf);


% 3. get the number of sampling points
max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2));


% 3. calculate sum of the curvature at sampling points
% Du                 1st-order derivative of the sampling points in u-direction
% Duu                2nd-order derivative of the sampling points in u-direction
% Dv                 1st-order derivative of the sampling points in v-direction
% Dvv                2nd-order derivative of the sampling points in v-direction
% ut_Fcur            discrete u in segments for the calculation of the curvature
% vt_Fcur            discrete v in segments for the calculation of the curvature

cur = zeros(1,numP_P);

if delta_v==0
   segType=1;
elseif delta_u==0
   segType=2;
end

if segType==1           % u-segment with constant v
    utFcur = linspace(TupleSeg(1,1),TupleSeg(2,1),numP_P);
    vtFcur = zeros(1,numP_P) + TupleSeg(1,2);
    
    [~, temp2, temp3] = nrbdeval(surf, surfD, surfDD, {utFcur vtFcur});
    temp4 = temp2{1,1};
    Du = temp4(:,:,1)';
    temp5 = temp3{1,1};
    Duu = temp5(:,:,1)';
        
    for j = 1:numP_P
        cur(1,j) = norm(cross(Du(j,:),Duu(j,:)))/norm(Du(j,:))^3;
    end
    
    sumCur = sum(cur);
    
elseif segType==2       % v-segment with constant u
    utFcur = zeros(1,numP_P) + TupleSeg(1,1);
    vtFcur = linspace(TupleSeg(1,2),TupleSeg(2,2),numP_P);

    [~, temp6, temp7] = nrbdeval(surf, surfD, surfDD, {utFcur vtFcur});
    temp8 = permute(temp6{1,2},[1 3 2]);
    Dv = 1000*temp8(:,:,1)';
    temp9 = permute(temp7{2,2},[1 3 2]);
    Dvv = 1000*temp9(:,:,1)';
    
    for j = 1:numP_P
        cur(1,j) = 1000*norm(cross(Dv(j,:),Dvv(j,:)))/norm(Dv(j,:))^3;
    end
    
    sumCur = sum(cur);
end
