r11_filted = xlsread('r11_filted.xlsx');
r11_npeaks =xlsread('least_squares.xlsx','n_peaks');
r11_para = xlsread('r11least_squares1.xlsx','para');

r11_g = zeros(800,871);
x=1:800;
for j=1:871
    if r11_npeaks(j) == 0
        for i=1:800
            r11_g(i,j)=NaN;
        end
    elseif r11_npeaks(j) == 1
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2);
    elseif r11_npeaks(j) == 2
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2) + r11_para(4,j).*exp(-((x-r11_para(5,j))./(r11_para(6,j).*sqrt(2))).^2);
    elseif r11_npeaks(j) == 3
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2) + r11_para(4,j).*exp(-((x-r11_para(5,j))./(r11_para(6,j).*sqrt(2))).^2) + r11_para(7,j).*exp(-((x-r11_para(8,j))./(r11_para(9,j).*sqrt(2))).^2);
    elseif r11_npeaks(j) == 4
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2) + r11_para(4,j).*exp(-((x-r11_para(5,j))./(r11_para(6,j).*sqrt(2))).^2) + r11_para(7,j).*exp(-((x-r11_para(8,j))./(r11_para(9,j).*sqrt(2))).^2) + r11_para(10,j).*exp(-((x-r11_para(11,j))./(r11_para(12,j).*sqrt(2))).^2);
    elseif r11_npeaks(j) == 5
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2) + r11_para(4,j).*exp(-((x-r11_para(5,j))./(r11_para(6,j).*sqrt(2))).^2) + r11_para(7,j).*exp(-((x-r11_para(8,j))./(r11_para(9,j).*sqrt(2))).^2) + r11_para(10,j).*exp(-((x-r11_para(11,j))./(r11_para(12,j).*sqrt(2))).^2) + r11_para(13,j).*exp(-((x-r11_para(14,j))./(r11_para(15,j).*sqrt(2))).^2);
    elseif r11_npeaks(j) == 6
        r11_g(:,j)= r11_para(1,j).*exp(-((x-r11_para(2,j))./(r11_para(3,j).*sqrt(2))).^2) + r11_para(4,j).*exp(-((x-r11_para(5,j))./(r11_para(6,j).*sqrt(2))).^2) + r11_para(7,j).*exp(-((x-r11_para(8,j))./(r11_para(9,j).*sqrt(2))).^2) + r11_para(10,j).*exp(-((x-r11_para(11,j))./(r11_para(12,j).*sqrt(2))).^2) + r11_para(13,j).*exp(-((x-r11_para(14,j))./(r11_para(15,j).*sqrt(2))).^2) + r11_para(16,j).*exp(-((x-r11_para(17,j))./(r11_para(18,j).*sqrt(2))).^2);
    end
end

xlswrite('r11G.xlsx',r11_g);

% Quality of fit
wdt = r11_filted - r11_g; 
mN=zeros(1,871);
sigmaN=zeros(1,871);
for j=1:871
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
sigmaWFi = zeros(1,871);
for j=1:871
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
        
    


    
        