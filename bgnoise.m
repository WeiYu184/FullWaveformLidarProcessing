% background noise
r11v_normal=xlsread('r11v_normal.xlsx');
%BgMean=xlsread('E:\wave33\prepare\BgMean.xlsx');
%BgSDEV=xlsread('E:\wave33\prepare\BgSDEV.xlsx');
%Tn=BgMean+4*BgSDEV;
r11v_Bg=zeros(800,871);
m150=zeros(1,871);
s150=zeros(1,871);
for j=1:871
    for i=1:150
        m150(j)=m150(j)+r11v_normal(i,j);
    end
     m150(j)= m150(j)/150;
    for i=1:150
        s150(j)=s150(j)+(r11v_normal(i,j)-m150(j))*(r11v_normal(i,j)-m150(j))/149;
    end
    s150(j)=sqrt(s150(j));
end

Tn=m150+4*s150;
for j=1:871
    for i=1:800
        if r11v_normal(i,j)<Tn(j)
            r11v_Bg(i,j)=0;
        else
            r11v_Bg(i,j)=r11v_normal(i,j)-Tn(j);
        end
    end
end
xlswrite('r11v_bgnoise.xlsx',r11v_Bg);