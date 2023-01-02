% Waveform normalisation
rec11_v=xlsread('E:\wave33\prepare\rec11_v.xlsx');

vsum=zeros(1,871); % 波形各采样点的电压值之和

for j=1:871 % for each waveform
    for i=1:800 % for each sample point
        if rec11_v(i,j)>0
            vsum(j)=vsum(j)+rec11_v(i,j);
        end
    end
end

r11v_normal=rec11_v./vsum;
for j=1:871
    if vsum(j)==0
        r11v_normal(:,j)=0;
    end
end

xlswrite('E:\wave33\prepare\r11v_normal_biggerthan0sum.xlsx',r11v_normal);