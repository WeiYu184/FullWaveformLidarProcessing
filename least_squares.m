% ��С���˷����
r11_filted = xlsread('r11_filted.xlsx');
in_peak = xlsread('findpeak.xlsx','n_peak');
il_peak = xlsread('findpeak.xlsx','l_peak');
ipeaks = xlsread('findpeak.xlsx','peaks');
isigmas = xlsread('findpeak.xlsx','sigmas');
% Tn = xlsread('E:\wave33\Gaussian\mN.xlsx');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% �������󣬴洢��������ֵ���м���
para = zeros(18,871);   % ����ֵ
ca_res = zeros(800,871); % ����ֵƫ���� y - f(p,x)
iterations = zeros(1,871); % ��������
r11_f = zeros(800,871); % �����Ĳ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TN = zeros(1,871);
Tn = zeros(1,871);
m150=zeros(1,871);
s150=zeros(1,871);
for j=1:871
    for i=1:150
        m150(j)=m150(j)+r11_filted(i,j);
    end
     m150(j)= m150(j)/150;
    for i=1:150
        s150(j)=s150(j)+(r11_filted(i,j)-m150(j))*(r11_filted(i,j)-m150(j))/149;
    end
    s150(j)=sqrt(s150(j));
    TN(j)=m150(j)+4*s150(j);
end


TN0=(TN~=0);           % ��0��λ��
TNnum0 = sum(TN0(:));  % ��0�ĸ���
TN(TN==0) = [];        % ȥ��0
TN=sum(TN)/TNnum0;
xlswrite('TN.xlsx',TN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��С�������
for j=1:871
%    Tn(j) = 1.5555e-07;
    Tn(j) = TN;
    if in_peak(:,j) == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif in_peak(:,j) == 6 % ����������
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j) ipeaks(5,j) il_peak(5,j) isigmas(5,j)  ipeaks(6,j) il_peak(6,j) isigmas(6,j)]; % λ�� ��� �����ʼֵ
        f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)) + p(13).* exp(-(x-p(14)).^2./(2.*p(15).^2)) + p(16).* exp(-(x-p(17)).^2./(2.*p(18).^2)));
        
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=p';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res;
        r11_f(:,j)= f(p,x);
        % ����������
        if abs(para(5,j)-para(2,j))<6 || abs(para(8,j)-para(5,j))<6 || abs(para(11,j)-para(8,j))<6 ||  abs(para(14,j)-para(11,j))<6 || abs(para(17,j)-para(14,j))<6 || abs(para(8,j)-para(2,j))<6 || abs(para(14,j)-para(8,j))<6 || abs(para(11,j)-para(5,j))<6
             para(:,j)=zeros(18,1);
             in_peak(:,j) = 5;
        end
        if para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2
             para(1,j)=0; para(2,j)=0; para(3,j)=0; 
        end
        if  para(4,j) < Tn(j) || para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2
            para(4,j)=0; para(5,j)=0; para(6,j)=0;
        end
        if para(7,j) < Tn(j) || para(8,j) < 0  ||  para(8,j) > 800 || para(9,j) < 2
            para(7,j)=0; para(8,j)=0; para(9,j)=0;
        end
        if para(10,j) < Tn(j) || para(11,j) < 0  ||  para(11,j) > 800 || para(12,j) < 2
            para(10,j)=0; para(11,j)=0; para(12,j)=0;
        end
        if para(13,j) < Tn(j) || para(14,j) < 0  ||  para(14,j) > 800 || para(15,j) < 2
            para(13,j)=0; para(14,j)=0; para(15,j)=0;
         end
        if para(16,j) < Tn(j) || para(17,j) < 0  ||  para(17,j) > 800 || para(18,j) < 2
            para(16,j)=0; para(17,j)=0; para(18,j)=0;
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    elseif in_peak(:,j) == 5 % ���������
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j) ipeaks(5,j) il_peak(5,j) isigmas(5,j)]; % λ�� ��� �����ʼֵ
        f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)) + p(13).* exp(-(x-p(14)).^2./(2.*p(15).^2)));
        
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=[p,0,0,0]';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res; 
        r11_f(:,j)= f(p,x);
        % ����������
        if abs(para(5,j)-para(2,j))<6 || abs(para(8,j)-para(5,j))<6 || abs(para(11,j)-para(8,j))<6 ||  abs(para(14,j)-para(11,j))<6 || abs(para(8,j)-para(2,j))<6 || abs(para(14,j)-para(8,j))<6 || abs(para(11,j)-para(5,j))<6
             para(:,j)=zeros(18,1);
             in_peak(:,j) = 4;
        end
        if para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2
             para(1,j)=0; para(2,j)=0; para(3,j)=0; 
        end
        if  para(4,j) < Tn(j) || para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2
            para(4,j)=0; para(5,j)=0; para(6,j)=0;
        end
        if para(7,j) < Tn(j) || para(8,j) < 0 ||  para(8,j) > 800 || para(9,j) < 2
            para(7,j)=0; para(8,j)=0; para(9,j)=0;
        end
        if para(10,j) < Tn(j) || para(11,j) < 0 ||  para(11,j) > 800 || para(12,j) < 2
            para(10,j)=0; para(11,j)=0; para(12,j)=0;
        end
         if para(13,j) < Tn(j) || para(14,j) < 0 ||  para(14,j) > 800 || para(15,j) < 2
            para(13,j)=0; para(14,j)=0; para(15,j)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif in_peak(:,j)== 4 % ���ĸ�����
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j)]; % λ�� ��� �����ʼֵ
        f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)));
        
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=[p,0,0,0,0,0,0]';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res;   
        r11_f(:,j)= f(p,x);
        % ����������
        if abs(para(5,j)-para(2,j))<6 || abs(para(8,j)-para(5,j))<6 || abs(para(11,j)-para(8,j))<6 || abs(para(8,j)-para(2,j))<6 || abs(para(11,j)-para(5,j))<6 || abs(para(11,j)-para(2,j))<6
             para(:,j)=zeros(18,1);
             in_peak(:,j) = 3;
        end
        if para(1,j) < Tn(j)|| para(2,j) < 0  ||  para(2,j) > 800 || para(3,j) < 2
             para(1,j)=0; para(2,j)=0; para(3,j)=0; 
        end
        if  para(4,j) < Tn(j) || para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2
            para(4,j)=0; para(5,j)=0; para(6,j)=0;
        end
        if para(7,j) < Tn(j) ||  para(8,j) < 0 ||  para(8,j) > 800 || para(9,j) < 2
            para(7,j)=0; para(8,j)=0; para(9,j)=0;
        end
        if para(10,j) < Tn(j) || para(11,j) < 0 ||  para(11,j) > 800 || para(12,j) < 2
            para(10,j)=0; para(11,j)=0; para(12,j)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif in_peak(:,j) == 3 % ����������
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j)]; % λ�� ��� �����ʼֵ
        f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)));
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=[p,0,0,0,0,0,0,0,0,0]';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res; 
        r11_f(:,j)= f(p,x);
        % ����������
        if abs(para(5,j)-para(2,j))<6 || abs(para(8,j)-para(5,j))<6 || abs(para(8,j)-para(2,j))<6
             para(:,j)=zeros(18,1);
             in_peak(:,j) = 2;
        end
        if para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2
             para(1,j)=0; para(2,j)=0; para(3,j)=0; 
        end
        if  para(4,j) < Tn(j) || para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2
            para(4,j)=0; para(5,j)=0; para(6,j)=0;
        end
        if para(7,j) < Tn(j) || para(8,j) < 0 ||  para(8,j) > 800 || para(9,j) < 2
            para(7,j)=0; para(8,j)=0; para(9,j)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif in_peak(:,j) == 2 % ����������
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j)]; % ��� λ�� �����ʼֵ
        f=@(p,x)(p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2)));
        
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=[p,0,0,0,0,0,0,0,0,0,0,0,0]';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res;  
        r11_f(:,j)= f(p,x);
        % 
        if abs(para(5,j)-para(2,j))<6
             para(:,j)=zeros(18,1);
             in_peak(:,j) = 1;
        end
        if  (para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2) && (para(4,j) < Tn(j) || para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2)
             para(:,j)=zeros(18,1);
        elseif para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2
             para(1,j)=0; para(2,j)=0; para(3,j)=0; 
        elseif  para(4,j) < Tn(j)|| para(5,j) < 0 ||  para(5,j) > 800 || para(6,j) < 2
            para(4,j)=0; para(5,j)=0; para(6,j)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif in_peak(:,j) == 1 % ��һ������
        y = r11_filted(:,j);
        x=1:800;x=x';
        p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j)]; % ��� λ�� ����ĳ�ʼֵ
        f=@(p,x)p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2));
        [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
        para(:,j)=[p,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
        iterations(j)=output.iterations;
        res = y - f(p,x);
        ca_res(:,j)=res;
        r11_f(:,j)= f(p,x);
        % ��������1
        if  para(1,j) < Tn(j) || para(2,j) < 0 ||  para(2,j) > 800 || para(3,j) < 2
            para(:,j)=zeros(18,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %elseif in_peak(:,j) == 2 % ����������
    %    y = r11_filted(:,j);
    %    x=1:800;x=x';
    %   p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j)]; % ��� λ�� �����ʼֵ
%         f=@(p,x)(p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2)));
%         
%         [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
%         para(:,j)=[p,0,0,0,0,0,0,0,0,0,0,0,0]';
%         iterations(j)=output.iterations;
%         res = y - f(p,x);
%         ca_res(:,j)=res; 
%          if  para(1,j) < Tn(j) || para(3,j) < 2
%             para(1,j)=0; para(2,j)=0; para(3,j)=0;
%          end
%          if  para(4,j) < Tn(j) || para(6,j) < 2
%             para(4,j)=0; para(5,j)=0; para(6,j)=0;
%          end
%          if abs(para(5,j)-para(2,j))<6
%              para(:,j)=zeros(18,1);
%              in_peak(:,j) = 1;
%          end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif in_peak(:,j) == 3 % ����������
%         y = r11_filted(:,j);
%         x=1:800;x=x';
%         p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j)]; % λ�� ��� �����ʼֵ
%         f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)));
%         [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
%         para(:,j)=[p,0,0,0,0,0,0,0,0,0]';
%         iterations(j)=output.iterations;
%         res = y - f(p,x);
%         ca_res(:,j)=res;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif in_peak(:,j)== 4 % ���ĸ�����
%         y = r11_filted(:,j);
%         x=1:800;x=x';
%         p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j)]; % λ�� ��� �����ʼֵ
%         f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)));
%         
%         [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
%         para(:,j)=[p,0,0,0,0,0,0]';
%         iterations(j)=output.iterations;
%         res = y - f(p,x);
%         ca_res(:,j)=res;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif in_peak(:,j) == 5 % ���������
%         y = r11_filted(:,j);
%         x=1:800;x=x';
%         p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j) ipeaks(5,j) il_peak(5,j) isigmas(5,j)]; % λ�� ��� �����ʼֵ
%         f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)) + p(13).* exp(-(x-p(14)).^2./(2.*p(15).^2)));
%         
%         [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
%         para(:,j)=[p,0,0,0]';
%         iterations(j)=output.iterations;
%         res = y - f(p,x);
%         ca_res(:,j)=res;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif in_peak(:,j) == 6 % ����������
%         y = r11_filted(:,j);
%         x=1:800;x=x';
%         p0=[ipeaks(1,j) il_peak(1,j) isigmas(1,j) ipeaks(2,j) il_peak(2,j) isigmas(2,j) ipeaks(3,j) il_peak(3,j) isigmas(3,j) ipeaks(4,j) il_peak(4,j) isigmas(4,j) ipeaks(5,j) il_peak(5,j) isigmas(5,j)  ipeaks(6,j) il_peak(6,j) isigmas(6,j)]; % λ�� ��� �����ʼֵ
%         f=@(p,x)( p(1).* exp(-(x-p(2)).^2./(2.*p(3).^2))+p(4).* exp(-(x-p(5)).^2./(2.*p(6).^2))+p(7).* exp(-(x-p(8)).^2./(2.*p(9).^2)) + p(10).* exp(-(x-p(11)).^2./(2.*p(12).^2)) + p(13).* exp(-(x-p(14)).^2./(2.*p(15).^2)) + p(16).* exp(-(x-p(17)).^2./(2.*p(18).^2)));
%         
%         [p,~,~,output] = fminsearch(@(p)sum((y-f(p,x)).*(y-f(p,x))),p0,optimset('MaxIter',10000000,'MaxFunEvals',10000000));
%         para(:,j)=[p]';
%         iterations(j)=output.iterations;
%         res = y - f(p,x);
%         ca_res(:,j)=res;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end

n_peaks = zeros(1,187);
for j=1:871
    no0=(para(:,j)~=0);   % ��0��λ��
    paranum = sum(no0(:));  % ��0�ĸ���
    n_peaks(j) = paranum/3;
end
xlswrite('least_squares.xlsx',para,'para');
xlswrite('least_squares.xlsx',ca_res,'ca_res');
xlswrite('least_squares.xlsx',iterations,'iterations');
xlswrite('least_squares.xlsx',n_peaks,'n_peaks');
xlswrite('r11_f.xlsx',r11_f);