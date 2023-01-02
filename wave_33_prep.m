% 波形预处理：量化到0-255 800个采样点*871
% piexl为各个坐标点的数据值（将其转换到0到255之间）
rec11=xlsread('rec11.xlsx');
rmax=255;% 要归一的范围的最大值
rmin=0;% 要归一的范围的最小值
r11max=max(max(rec11));% 所有数据中最大的
r11min=min(min(rec11));% 所有数据中最小的
rec11_255=zeros(800,871);

for j=1:871
    for i=1:800
        rec11_255(i,j) = round((rmax-rmin)*(rec11(i,j)-r11min)/(r11max-r11min) + rmin); %归一化并取整
    end
end

xlswrite('rec11_255.xlsx',rec11_255);

% look-up table.xlsx电压值转换
vtable=xlsread('E:\wave33\prepare\look-up table.xlsx');
rec11_v=zeros(800,871);
for j=1:871
    for i=1:800
        v_count=rec11_255(i,j)+1;
        rec11_v(i,j)=vtable(v_count);
    end
end
xlswrite('rec11_v.xlsx',rec11_v);