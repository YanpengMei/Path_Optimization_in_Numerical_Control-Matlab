% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% combine segments to get the path

function Path = A5_Seg2Path(Seg)

for i= 1:size(Seg,2)
    temp1 = Seg{i};
    if i==1
        Path = temp1;
    else
        temp2 = temp1(:,2:end);
        Path = [Path temp2];
    end
end

end