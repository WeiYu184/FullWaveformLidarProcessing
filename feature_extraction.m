% ����������ȡ
para=xlsread('least_squares.xlsx','para');
n_peaks=xlsread('least_squares.xlsx','n_peaks');
%r11_g = xlsread('r11_f.xlsx');
mN = xlsread('mN.xlsx');
TN = xlsread('TN.xlsx');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����С����������ݽ��в�0
% ȥ���м��0����ĩβ��0
for j=1:871
    a= para(:,j);
    a0=(a==0);           % 0��λ��
    anum0 = sum(a0(:));  % 0�ĸ���
    a(a==0) = [];        % ȥ��0
    para(:,j) =[a;zeros(anum0,1)]; % ĩβ��0
end
xlswrite('least_squares1.xlsx',para,'para');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���㲨����������
pEcho = zeros(1,871);    % ����ȫ�ߣ���Ч�źŷ�Χ����_�ڲ�߶�
wpeak = zeros(1,871);    % ���γ��ȣ��źſ�ʼ��������һ����ֵ�ĳ���_�ڲ�߶�
dpeak = zeros(1,871);    % ��ֵ���룺��һ����ֵ��������һ����ֵ�ĳ���_�ڲ�߶�

pbeg = zeros(1,871); % �������
pend = zeros(1,871); % �����յ�
peak0 = zeros(1,871);% ��һ�������λ��
peakG = zeros(1,871);% ���һ�������λ��

pCenter = zeros(1,871);% ��������
HOME = zeros(1,871);% �벨������

for j=1:871
    if n_peaks(j)~= 0
        l_peakG = 3*n_peaks(j)-1;
        peakG(j) = para(l_peakG,j);
        peak0(j) = para(2,j);
    end
    pbeg(j) = 1;
    pend(j) = 800;
    noBg=(r11_g(:,j)>TN);   % ��Ч�źŵ�λ��
    % ����pbeg��pend��
    while noBg(pbeg(j)) == 0 && pbeg(j) < 800
        pbeg(j)=pbeg(j)+1;
    end
    while noBg(pend(j)) == 0 && pend(j) > 1
        pend(j)=pend(j)-1;
    end
    % ���㲨��ȫ�ߣ�
    pEcho(j)=pend(j) - pbeg(j);
    if abs(pEcho(j))> 798
        pEcho(j)= 0;pend(j)=0; pbeg(j)=0;
    end
    % �����ֵ���룺
    dpeak(j)= abs(peakG(j)- peak0(j));
    % ���㲨�γ���:
    if peakG(j)> pbeg(j)
        wpeak(j) = peakG(j) - pbeg(j);
    end
    % ����벨�����ߣ�
    ensum = sum(r11_g(:,j));
    ensum50 = ensum/2;
    en_sum = 0;
    i = 0;
    while en_sum < ensum50
        i = i + 1;
        en_sum = en_sum + r11_g(i,j);
    end
    pCenter(j) = i;
    HOME(j) = peakG(j)-pCenter(j);
end  
    
xlswrite('feature_para.xlsx',pEcho,'����ȫ��');
xlswrite('feature_para.xlsx',wpeak,'wpeak');
xlswrite('feature_para.xlsx',dpeak,'dpeak');
xlswrite('feature_para.xlsx',pbeg,'�������');
xlswrite('feature_para.xlsx',pend,'�����յ�');
xlswrite('feature_para.xlsx',peakG,'peakG');
xlswrite('feature_para.xlsx',peak0,'peak0');
xlswrite('feature_para.xlsx',pCenter,'pCenter');
xlswrite('feature_para.xlsx',HOME,'HOME');
xlswrite('feature_para.xlsx',n_peaks,'������');