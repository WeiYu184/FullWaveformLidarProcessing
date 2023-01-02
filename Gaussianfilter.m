% 功能：对一维信号的高斯滤波，头尾r/2的信号不进行滤波
% r     :高斯模板的大小推荐奇数
% sigma :标准差
% y     :需要进行高斯滤波的序列
%a1 =     0.01712  (0.01702, 0.01722) 0.06275  (0.06238, 0.06312)
%b1 =       217.8  (217.8, 217.9)
%c1 =        9.12  (9.058, 9.182)
 
sigma=9.12/sqrt(2);
r=218;
tmp11_normal=xlsread('tmp11_normal.xlsx');
r11v_Bg=xlsread('r11v_bgnoise.xlsx');
% 生成一维高斯滤波模板
GaussTemp11 = zeros(1,r*2-1);
for i=1 : r*2-1
   GaussTemp11(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end

% 高斯滤波 r11v_Bg → r11_filt:
r11_filt=zeros(1234,871);
for j=1:871
    r11_filt(:,j) = conv(r11v_Bg(:,j),GaussTemp11);
end
r11_filted=r11_filt(218:1017,:);
xlswrite('r11_filted.xlsx',r11_filted);
%for i = r : length(y)-r
%   y_filted(i) = y(i-r+1 : i+r-1)*GaussTemp11';
%end