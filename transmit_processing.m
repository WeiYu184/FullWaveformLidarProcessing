% transmit11 processing
tra11=xlsread('tra11.xlsx');
% ����Ԥ����������0-255 800��������*871
% piexlΪ��������������ֵ������ת����0��255֮�䣩
tmax=255;% Ҫ��һ�ķ�Χ�����ֵ
tmin=0;% Ҫ��һ�ķ�Χ����Сֵ
t11max=1023;% ��������������
t11min=0;% ������������С��
tra11_255=zeros(400,871);

for j=1:871
    for i=1:400
        tra11_255(i,j) = round((tmax-tmin)*(tra11(i,j)-t11min)/(t11max-t11min) + tmin); %��һ����ȡ��
    end
end

% look-up table.xlsx��ѹֵת��
vtable=xlsread('E:\wave33\prepare\look-up table.xlsx');
tra11_v=zeros(400,871);
for j=1:871
    for i=1:400
        v_count=tra11_255(i,j)+1;
        tra11_v(i,j)=vtable(v_count);
    end
end

% Waveform normalisation
vsum=sum(tra11_v); % ���θ�������ĵ�ѹֵ֮��
t11v_normal=zeros(400,871); % ���ι�һ����ĵ�ѹֵ

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

% ȡƽ�����õ���˹�˲�ģ��
tmp11 = mean(t11v_Bg,2);
xlswrite('E:\wave33-11\tmp11.xlsx',tmp11);

% Waveform normalisation
tmp11sum=sum(tmp11); % ���θ�������ĵ�ѹֵ֮��
tmp11_normal=tmp11./tmp11sum;

xlswrite('E:\wave33-11\tmp11_normal.xlsx',tmp11_normal);