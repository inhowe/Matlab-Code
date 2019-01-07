% 长除法 
% 参数分别为分子向量、分母向量、结果的长度（默认为5）

function res = longDiv(nom, den, bit)
    if nargin < 3
        bit = length(den) * 2;
    end
    if length(den) < length(nom)
        disp('error z transform');
        return;
    end
    if length(den) ~= length(nom)
       nom = [zeros(1, length(den) - length(nom))   , nom]; 
    end

    res = [];
    m = nom;
    for i = 1 : bit
        tempRes = m(1)/den(1);
        m = m - tempRes * den;
        m = [m(2:length(m)), 0];
        res = [res tempRes];
    end
end
