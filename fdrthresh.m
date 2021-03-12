function thresh = fdrthresh(z,q)

az = abs(z);
[sz,iz] = sort(az);
pobs = erfc(sz./sqrt(2));
N = 1:length(z);
pnull =  N' ./length(z);
good = (fliplr(pobs) <= (q .* pnull));
if any(good),
    izmax  = iz(1 + length(sz) - max(N(good)));
    thresh = az(izmax);
else
    thresh = max(az) + .01;
end
    
%
% Copyright (c) 2006. Iddo Drori
%  

