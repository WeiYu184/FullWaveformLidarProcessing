% Quality of fit
r11_g = xlsread('r11G.xlsx');
r11_F = xlsread('r11_filted.xlsx');
wdt = r11_F - r11_g; 
mN=zeros(1,845);
sigmaN=zeros(1,845);
for j=1:845
    for i=1:150
        mN(j)=mN(j)+r11_g(i,j);
    end
     mN(j)= mN(j)/150;
    for i=1:150
        sigmaN(j)=sigmaN(j)+(r11_g(i,j)-mN(j))*(r11_g(i,j)-mN(j));
    end
    sigmaN(j)=sqrt(sigmaN(j)/149); % sigmaN
end
sigma_N = sum(sigmaN);

mWFi = sum(wdt);
mWFi = mWFi./800;
sigmaWFi = zeros(1,845);
for j=1:845
    for i=1:800
    sigmaWFi(j) = sigmaWFi(j) + (wdt(i,j)-mWFi(j))*(wdt(i,j)-mWFi(j));
    end
    sigmaWFi(j) = sqrt(sigmaWFi(j)/799);
end

sigma_WF = sum(sigmaWFi)/871;
sigmaWFi_sigmaN = sigmaWFi./sigmaN;
sigmaWFi_sigma_N = sigmaWFi./sigma_N;

xlswrite('r11sigmaWF.xlsx', sigmaN,'sigmaN');
xlswrite('r11sigmaWF.xlsx',sigmaWFi,'sigmaWFi');