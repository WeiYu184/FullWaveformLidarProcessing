% transmit11 processing
tra11=xlsread('tra11.xlsx');
% 波形预处理：量化到0-255 800个采样点*871
% piexl为各个坐标点的数据值（将其转换到0到255之间）
tmax=255;% 要归一的范围的最大值
tmin=0;% 要归一的范围的最小值
t11max=1023;% 所有数据中最大的
t11min=0;% 所有数据中最小的
tra11_255=zeros(400,871);

for j=1:871
    for i=1:400
        tra11_255(i,j) = round((tmax-tmin)*(tra11(i,j)-t11min)/(t11max-t11min) + tmin); %归一化并取整
    end
end

% look-up table.xlsx电压值转换
vtable=xlsread('E:\wave33\prepare\look-up table.xlsx');
tra11_v=zeros(400,871);
for j=1:871
    for i=1:400
        v_count=tra11_255(i,j)+1;
        tra11_v(i,j)=vtable(v_count);
    end
end

% Waveform normalisation
vsum=sum(tra11_v); % 波形各采样点的电压值之和
t11v_normal=zeros(400,871); % 波形归一化后的电压值

for j=1:871 % for each waveform
    for i=1:400 % for each sample point
        t11v_normal(i,j)=tra11_v(i,j)/abs(vsum(j));
    end
end

%xlswrite('E:\wave33-11\t11v_normal.xlsx',t11v_normal);

% background noise
t11v_Bg=zeros(400,871);
m150=zeros(1,871);
s150=zeros(1,871);
for j=1:871
    for i=1:150
        m150(j)=m150(j)+t11v_normal(i,j);
    end
     m150(j)= m150(j)/150;
    for i=1:150
        s150(j)=s150(j)+(t11v_normal(i,j)-m150(j))*(t11v_normal(i,j)-m150(j))/149;
    end
    s150(j)=sqrt(s150(j));
end

Tn=m150+4*s150;
for j=1:871
    for i=1:400
        if t11v_normal(i,j)<Tn(j)
            t11v_Bg(i,j)=0;
        else
            t11v_Bg(i,j)=t11v_normal(i,j)-Tn(j);
        end
    end
end
%xlswrite('E:\wave33-11\t11v_bgnoise.xlsx',t11v_Bg);

% 取平均，得到高斯滤波模板
tmp11 = mean(t11v_Bg,2);
xlswrite('E:\wave33-11\tmp11.xlsx',tmp11);

% Waveform normalisation
tmp11sum=sum(tmp11); % 波形各采样点的电压值之和
tmp11_normal=tmp11./tmp11sum;

xlswrite('E:\wave33-11\tmp11_normal.xlsx',tmp11_normal);