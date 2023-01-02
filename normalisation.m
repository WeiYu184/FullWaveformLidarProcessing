% Waveform normalisation
rec11_v=xlsread('rec11_v.xlsx');
vsum=sum(rec11_v); % 波形各采样点的电压值之和
r11v_normal=zeros(n_v,n_waveform); % 波形归一化后的电压值

for j=1:n_waveform % for each waveform
    for i=1:n_v % for each sample point
        r11v_normal(i,j)=rec11_v(i,j)/abs(vsum(j));
    end
end
xlswrite('r11v_normal.xlsx',r11v_normal);
