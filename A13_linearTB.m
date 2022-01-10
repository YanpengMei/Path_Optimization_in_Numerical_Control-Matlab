% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% define the linear tolerance band
% input of this function should be the control points of the real surface !!!


function [disLTB2RP] = A13_linearTB(cp_var)
load('A2_parameters.mat','tb_r','knot')
load('A2_fittedCP_best.mat','fittedCP_best')


% 1. get the real path in 2D
[~,~,xx,yy] = A6_fit(cp_var);
realpath2D = [xx ; yy];


% 2. define and explain the parameters
% numPP              number of pieces of the original path
numPP = size(realpath2D,2) - 1;

% numPb2P           number of points between each 2 points of the path
% note that, this "numPb2P" includes the points of the real path
numPb2P = 10;

% numP_fine         number of points in finer path
numP_fine = numPb2P*numPP-(numPP-1);


% 3. get finer real path in 2D(real path with more sampling points)
% path_fine         real path with more sampling points
curve = nrbmak(fittedCP_best,knot);
ut_path_fine = linspace(0,1,numP_fine);
temp1 = nrbeval(curve,ut_path_fine);
temp2 = temp1(:,:,1);
path_fine = temp2(1:2,:);


% 4. get the pieces of the finer path
% 'finer' means more sampling points
path_pieces = zeros(2,numPb2P,numPP);
for i = 1:numPP
    if i ==1
        path_pieces(:,:,i) = path_fine(:,1:numPb2P);
    else
        temp3 = numPb2P*(i-1)-(i-2);
        temp4 = temp3 + (numPb2P - 1);
        path_pieces(:,:,i) = path_fine(:,temp3:temp4);
    end
end


% 5. build the line-function "y=k*x+b" for each piece
delta_y = realpath2D(2,2:end) - realpath2D(2,1:(end-1));
delta_x = realpath2D(1,2:end) - realpath2D(1,1:(end-1));
k = delta_y./delta_x;
% from "y0=k*x0+b" we know that "b=y0-k*x0"
x0 = realpath2D(1,1:(end-1));
y0 = realpath2D(2,1:(end-1));
b = y0-k.*x0;


% 6. calculate the minimum distance between the real points and the line-function for each piece
% if we write the line-function in the form Ax+By+C=0, then,
A = k;
B = zeros(1,numPP) - 1;
C = b;

% disPL            distance between the real points and the line-function
% disPL_max        maximum distance between the real points and the line-function
disPL = zeros(1,numPb2P);
disPL_max = zeros(1,numPP);

for i=1:numPP
    for j=1:numPb2P
        disPL(j) = abs( A(i)*path_pieces(1,j,i) + B(i)*path_pieces(2,j,i) + C(i))/sqrt(A(i)^2 + B(i)^2);
    end
    disPL_max(i) = max(disPL);
end


% 7. get parameters for tb_up and tb_down
% disLTB2RP         distance between linear tolerance band and real path
disLTB2RP = tb_r - disPL_max;



end