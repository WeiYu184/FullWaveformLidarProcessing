% ����Ԥ����������0-255 800��������*871
% piexlΪ��������������ֵ������ת����0��255֮�䣩
rec11=xlsread('rec11.xlsx');
rmax=255;% Ҫ��һ�ķ�Χ�����ֵ
rmin=0;% Ҫ��һ�ķ�Χ����Сֵ
r11max=max(max(rec11));% ��������������
r11min=min(min(rec11));% ������������С��
rec11_255=zeros(800,871);

for j=1:871
    for i=1:800
        rec11_255(i,j) = round((rmax-rmin)*(rec11(i,j)-r11min)/(r11max-r11min) + rmin); %��һ����ȡ��
    end
end

xlswrite('rec11_255.xlsx',rec11_255);

% look-up table.xlsx��ѹֵת��
vtable=xlsread('E:\wave33\prepare\look-up table.xlsx');
rec11_v=zeros(800,871);
for j=1:871
    for i=1:800
        v_count=rec11_255(i,j)+1;
        rec11_v(i,j)=vtable(v_count);
    end
end
xlswrite('rec11_v.xlsx',rec11_v);