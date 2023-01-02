% 初始参数估计
r11_filted=xlsread('r11_filted.xlsx');
n_peak=zeros(1,871);   %峰值数量
l_peak=zeros(6,871); %峰值位置
peaks=zeros(6,871);  %峰值大小
sigmas=zeros(6,871);  %脉宽
for j=1:871 % for each waveform
    for i=1:800-4 % 遍历搜索
        win=[r11_filted(i,j) r11_filted(i+1,j) r11_filted(i+2,j) r11_filted(i+3,j) r11_filted(i+4,j)]; %搜索窗
        if (max(win) == r11_filted(i+2,j)) && (min(win) < r11_filted(i+2,j))
            n_peak(1,j)=n_peak(1,j)+1;
            n=n_peak(1,j);
            peaks(n,j)= r11_filted(i+2,j);
            l_peak(n,j)=i+2;
        elseif (min(r11_filted(i+1,j),r11_filted(i+3,j)) > max(r11_filted(i,j),r11_filted(i+4,j)))
            n_peak(1,j)=n_peak(1,j)+1;
            n=n_peak(1,j);
            peaks(n,j)= r11_filted(i+2,j);
            l_peak(n,j)=i+2;
        end
     end
end

% 合并间隔小于10的峰值
for j=1:871
    for i=1:5
        if (l_peak(i+1,j)-l_peak(i,j))<10 && (l_peak(i+1,j)-l_peak(i,j))>0
            n_peak(1,j)=n_peak(1,j)-1;
            l_peak(i+1,j)=round((l_peak(i,j)+l_peak(i+1,j))/2);
            l_peak(i,j)=0;
            peaks(i+1,j)=max(peaks(i,j),peaks(i+1,j));
            peaks(i,j)=0;          
        end        
    end
    % 去除中间的0并在末尾补0
    a= peaks(:,j);
    a0=(a==0);           % 0的位置
    anum0 = sum(a0(:));  % 0的个数
    a(a==0) = [];        % 去掉0
    peaks(:,j) =[a;zeros(anum0,1)]; % 末尾补0
    b= l_peak(:,j);
    b0=(b==0);           % 0的位置
    bnum0 = sum(b0(:));  % 0的个数
    b(b==0) = [];        % 去掉0
    l_peak(:,j) =[b;zeros(bnum0,1)]; % 末尾补0
end

% 求脉宽
pmax = max(peaks);
for j=1:871
    if n_peak(1,j) == 0
        sigmas(1,j)=0;
    elseif n_peak(1,j) == 1
        sigmas(1,j)=9.121/sqrt(2);
    elseif peaks(1,j) == pmax(j)
         sigmas(1,j)=9.121/sqrt(2);
    else
        sigmas(1,j)= (l_peak(2,j)-l_peak(1,j))/2;
    end  
    for i=2:6
        if l_peak(i,j) ~=0 &&  peaks(i,j) == pmax(j)
            sigmas(i,j)=9.121/sqrt(2);
        elseif l_peak(i,j) ~=0
            sigmas(i,j)=(l_peak(i,j)-l_peak(i-1,j))/2;
        end
    end
end
xlswrite('findpeak.xlsx',n_peak,'n_peak');
xlswrite('findpeak.xlsx',l_peak,'l_peak');
xlswrite('findpeak.xlsx',peaks,'peaks');
xlswrite('findpeak.xlsx',sigmas,'sigmas');

%绘图
x = 1:1:800;
y = exp(-(x-184).*(x-184)/(2.*15.5.*15.5)).*0.0029;
plot(x,y);

x = 1:1:800;
y = exp(-(x-215).*(x-215)/(2.*15.5.*15.5)).*0.0063;
plot(x,y);

x = 1:1:800;
y = exp(-(x-250).*(x-250)/(2.*17.5.*17.5)).*0.0072;
plot(x,y);

x = 1:1:800;
y = exp(-(x-361).*(x-361)/(2.*6.4495.*6.4495)).*0.0249;
plot(x,y);
