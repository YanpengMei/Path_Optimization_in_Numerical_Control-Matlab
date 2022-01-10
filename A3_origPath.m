% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% Get a discreted path with given Tuple, sampling density and test surface
% Note: this path is the original path, input should be 'cp'.


function path = A3_origPath(cp_var)

load('A2_parameters.mat','knots','Tuple','cp','step_xy','segType','numSeg')
surface = nrbmak(cp_var,knots);

% get the "range" of the surface by measuring the control points
% max_delta_cpX        maximum difference of x-coordinate of each row of control points
% max_delta_cpY        maximum difference of y-coordinate of each column of control points
max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));

% dis_Tuple            distance of the points in Tuple
delta_u = (Tuple(2:end,1) - Tuple(1:(end-1),1));
delta_v = (Tuple(2:end,2) - Tuple(1:(end-1),2));

% numP_P               number of points in each segments of the path
numP_P = ceil(sqrt(((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2)));


% get the path
% path                  generated path
% numP_path             number of points in the path
numP_path = sum(numP_P) - (numSeg-1);
path = zeros(3,numP_path); 

% ut_path             discrete u for the path
% vt_path             discrete v for the path

temp1 = 0;

% we can devided the path generation into 3 cases:
% case 1: iso-parametric segments with constant v
% case 2: iso-parametric segments with constant u
% case 3: non-iso-parametric segments

for i = 1:numSeg
    if segType(i,2)==1       % u-segment with constant v
        ut_path = linspace(Tuple(i,1),Tuple(i+1,1),numP_P(i));
        vt_path = zeros(1,numP_P(i)) + Tuple(i,2);  
        temp2 = nrbeval(surface,{ut_path,vt_path});
        if i==1
            temp1 = temp1 + numP_P(i);
            path(:,1:numP_P(i)) = temp2(:,:,1);
        else
            temp1 = temp1 + numP_P(i) - 1;
            path(:,(temp1 - numP_P(i) + 1):temp1) = temp2(:,:,1);
        end
        
    elseif segType(i,2)==2       % v-segment with constant u
        ut_path = zeros(1,numP_P(i)) + Tuple(i,1);
        vt_path = linspace(Tuple(i,2),Tuple(i+1,2),numP_P(i));
        temp3 = nrbeval(surface,{ut_path,vt_path});
        temp4 = permute(temp3,[1 3 2]);
        if i==1
            temp1 = temp1 + numP_P(i);
            path(:,1:numP_P(i)) = temp4(:,:,1);
        else
            temp1 = temp1 + numP_P(i) - 1;
            path(:,(temp1 - numP_P(i) + 1):temp1) = temp4(:,:,1);
        end
        
    else % 
        ut_path = linspace(Tuple(i,1),Tuple(i+1,1),numP_P(i));
        vt_path = linspace(Tuple(i,2),Tuple(i+1,2),numP_P(i));
        temp5 = zeros(3,numP_P(i));
        for j = 1:numP_P(i)
            temp5(:,j) = nrbeval(surface,{ut_path(j),vt_path(j)});
        end
        if i==1
            temp1 = temp1 + numP_P(i);
            path(:,1:numP_P(i)) = temp5(:,:,1);
        else
            temp1 = temp1 + numP_P(i) - 1;
            path(:,(temp1 - numP_P(i) + 1):temp1) = temp5(:,:,1);
        end
    end
end
% end of the function main body
end