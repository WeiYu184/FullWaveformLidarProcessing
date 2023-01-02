% Waveform normalisation
rec11_v=xlsread('rec11_v.xlsx');
vsum=sum(rec11_v); % ���θ�������ĵ�ѹֵ֮��
r11v_normal=zeros(800,871); % ���ι�һ����ĵ�ѹֵ

for j=1:871 % for each waveform
    for i=1:800 % for each sample point
        r11v_normal(i,j)=rec11_v(i,j)/abs(vsum(j));
    end
end
xlswrite('r11v_normal.xlsx',r11v_normal);