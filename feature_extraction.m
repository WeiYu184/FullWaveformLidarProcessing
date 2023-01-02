% 特征参数提取
para=xlsread('least_squares.xlsx','para');
n_peaks=xlsread('least_squares.xlsx','n_peaks');
%r11_g = xlsread('r11_f.xlsx');
mN = xlsread('mN.xlsx');
TN = xlsread('TN.xlsx');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对最小二乘拟合数据进行补0
% 去除中间的0并在末尾补0
for j=1:871
    a= para(:,j);
    a0=(a==0);           % 0的位置
    anum0 = sum(a0(:));  % 0的个数
    a(a==0) = [];        % 去掉0
    para(:,j) =[a;zeros(anum0,1)]; % 末尾补0
end
xlswrite('least_squares1.xlsx',para,'para');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算波形特征参数
pEcho = zeros(1,871);    % 波形全高：有效信号范围长度_冠层高度
wpeak = zeros(1,871);    % 波形长度：信号开始到倒数第一个峰值的长度_冠层高度
dpeak = zeros(1,871);    % 峰值距离：第一个峰值到倒数第一个峰值的长度_冠层高度

pbeg = zeros(1,871); % 波形起点
pend = zeros(1,871); % 波形终点
peak0 = zeros(1,871);% 第一个波峰的位置
peakG = zeros(1,871);% 最后一个波峰的位置

pCenter = zeros(1,871);% 能量中心
HOME = zeros(1,871);% 半波能量高

for j=1:871
    if n_peaks(j)~= 0
        l_peakG = 3*n_peaks(j)-1;
        peakG(j) = para(l_peakG,j);
        peak0(j) = para(2,j);
    end
    pbeg(j) = 1;
    pend(j) = 800;
    noBg=(r11_g(:,j)>TN);   % 有效信号的位置
    % 计算pbeg和pend：
    while noBg(pbeg(j)) == 0 && pbeg(j) < 800
        pbeg(j)=pbeg(j)+1;
    end
    while noBg(pend(j)) == 0 && pend(j) > 1
        pend(j)=pend(j)-1;
    end
    % 计算波形全高：
    pEcho(j)=pend(j) - pbeg(j);
    if abs(pEcho(j))> 798
        pEcho(j)= 0;pend(j)=0; pbeg(j)=0;
    end
    % 计算峰值距离：
    dpeak(j)= abs(peakG(j)- peak0(j));
    % 计算波形长度:
    if peakG(j)> pbeg(j)
        wpeak(j) = peakG(j) - pbeg(j);
    end
    % 计算半波能量高：
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
    
xlswrite('feature_para.xlsx',pEcho,'波形全高');
xlswrite('feature_para.xlsx',wpeak,'wpeak');
xlswrite('feature_para.xlsx',dpeak,'dpeak');
xlswrite('feature_para.xlsx',pbeg,'波形起点');
xlswrite('feature_para.xlsx',pend,'波形终点');
xlswrite('feature_para.xlsx',peakG,'peakG');
xlswrite('feature_para.xlsx',peak0,'peak0');
xlswrite('feature_para.xlsx',pCenter,'pCenter');
xlswrite('feature_para.xlsx',HOME,'HOME');
xlswrite('feature_para.xlsx',n_peaks,'波峰数');