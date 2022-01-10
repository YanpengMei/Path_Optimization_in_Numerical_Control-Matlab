% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% to get the initial control points for fitting the discreted points into a NURBS-curve
% Note: first and last control points have also the information of coordinate transformation


function [CP_init,numCP,xx,yy] = A6_fit(cp_var)
load('A2_TupleSeg.mat','TupleSeg')


%% 1. get real segment in 2D
% 1. segment type
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));

% segType        segment type,"1" for constant v, "2" for constant u, "3" for non-iso...
if delta_v==0        % u-segment
        segType=1;
elseif delta_u==0    % v-segment
        segType=2;
else                 % non-iso
        segType=3;
end


% 2. get 2D real segment
if segType == 1
    temp1 = A4_getSeg(cp_var);
    xx = temp1(1,:) - temp1(1,1);
    yy = temp1(3,:) - temp1(3,1);
elseif segType == 2
    temp2 = A4_getSeg(cp_var);
    xx = temp2(2,:) - temp2(2,1);
    yy = temp2(3,:) - temp2(3,1);
else    
    temp3 = A4_getSeg(cp_var);
    xx = temp3(1,:);
    yy = temp3(2,:);
end



%% 2. get the initial fitted NURBS-curve                   
% 1. numCP             number of control points
% we use 14 control points to fit the NURBS-curve(14 is not optimal, but works well)
numCP = 14;


% 2. set initial control points
% first and last control points are also the first and last points of the segment
CP_init = zeros(2,numCP);
CP_init(1,1) = xx(1);
CP_init(2,1) = yy(1);
CP_init(1,numCP) = xx(end);
CP_init(2,numCP) = yy(end);

% the rest initial control points will be distributed uniformly
% delta_CPxx          difference of xx between first and last points of the segment
% delta_CPyy          difference of yy between first and last points of the segment
delta_CPxx = xx(end);
delta_CPyy = yy(end);

%const           a constant number to avoid linear distribution of control points
const = 0; % so far not used

for i = 2:numCP-1
    if rem(i,2)==0
        CP_init(1,i) = delta_CPxx*((i-1)/(numCP-1));
        CP_init(2,i) = delta_CPyy*((i-1)/(numCP-1))+const;
    else
        CP_init(1,i) = delta_CPxx*((i-1)/(numCP-1));
        CP_init(2,i) = delta_CPyy*((i-1)/(numCP-1))-const;
    end
end


end




