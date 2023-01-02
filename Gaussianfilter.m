% ���ܣ���һά�źŵĸ�˹�˲���ͷβr/2���źŲ������˲�
% r     :��˹ģ��Ĵ�С�Ƽ�����
% sigma :��׼��
% y     :��Ҫ���и�˹�˲�������
%a1 =     0.01712  (0.01702, 0.01722) 0.06275  (0.06238, 0.06312)
%b1 =       217.8  (217.8, 217.9)
%c1 =        9.12  (9.058, 9.182)
 
sigma=9.12/sqrt(2);
r=218;
tmp11_normal=xlsread('tmp11_normal.xlsx');
r11v_Bg=xlsread('r11v_bgnoise.xlsx');
% ����һά��˹�˲�ģ��
GaussTemp11 = zeros(1,r*2-1);
for i=1 : r*2-1
   GaussTemp11(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end

% ��˹�˲� r11v_Bg �� r11_filt:
r11_filt=zeros(1234,871);
for j=1:871
    r11_filt(:,j) = conv(r11v_Bg(:,j),GaussTemp11);
end
r11_filted=r11_filt(218:1017,:);
xlswrite('r11_filted.xlsx',r11_filted);
%for i = r : length(y)-r
%   y_filted(i) = y(i-r+1 : i+r-1)*GaussTemp11';
%end