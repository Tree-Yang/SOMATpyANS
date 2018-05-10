clear;clc;close all;
[history,timeconsumption,output] = runfmincon;
%clc;    %���������������Ϣ
disp('=========================================');
disp('             ����ִ����Ϣ��');
disp('-----------------------------------------')
disp(['�ܵ���������           ',num2str(output.iterations)]);
disp(['Ŀ�꺯��ִ�д�����     ',num2str(output.funcCount)]);
disp(['����������ʱ�䣺       ',num2str(timeconsumption),'  s']);
disp('���ŵ㣺');
disp('++++++++++++++++++++++');
disp(history.x(:,end));
disp('++++++++++++++++++++++');
disp(['���ź���ֵ��            ',num2str(history.fval(end))]);
disp('=========================================')
figure(1);
plot(history.fval,'-o','LineWidth',1.5);
grid on;
xlabel('Iteration');
ylabel('Objective Function');
set(gca,'FontName','Euclid','FontSize',18);
figure(2);
plot(history.x(1,:),'-o','LineWidth',1.5);
hold on;
plot(history.x(2,:),'--s','LineWidth',1.5);
plot(history.x(3,:),'-.d','LineWidth',1.5);
legend('\itA_{\rm1}','\itA_{\rm2}','\itA_{\rm3}','Location','SouthWest');
grid on;
xlabel('Iteration');
ylabel('Section Area');
set(gca,'FontName','Euclid','FontSize',18);
%===========================================================================================
function [history,timeconsumption,output] = runfmincon
%����һ���ְ����ں�����ֻ��Ϊ�����fmincon�����Ľ������ŵ�����    
    x0 = [0.5;0.5;0.5];%��ʼ��
    lb = [0.001;0.001;0.001];
    ub = [Inf;Inf;Inf];
    history.x = [];
    %+++++++++++++++++++++
    %��ͨ��output function������ź���ֵ���л�������������ʱ����Ҫ��������������ȡ��ע��
    history.fval = [];
    %searchdir=[];
    %++++++++++++++++����ݶ�++++++++++++++++
    options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed');  %fmincon����������
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed');  %fmincon����������
    %++++++++++++++++�����ݶ�++++++++++++++++
    %sqp�㷨
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed',...
    %    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);  %fmincon����������
    %���������㷨
    % options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter-detailed',...
    %     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);  %fmincon����������
    %++++++++++++++++++++++++++++++++++++
    tstart = tic;   %��ʼ��ʱ
    %[x,fval,exitflag,output,lambda] = fmincon(@pnonpobj,x0,[],[],[],[],lb,ub,@pnonpcon,options);
    [~,~,exitflag,output,~] = fmincon(@pnonpobj,x0,[],[],[],[],lb,ub,@pnonpcon,options);
    if exitflag ~= 1
        disp(['exitflag =  ',num2str(exitflag)]);
        error('fmincon����,����exitflag��ֵ����ԭ��!');
    end
    tend = toc(tstart);     %��ʱ����
    timeconsumption = tend;
%========================================================================================
    %�ֲ�����
%-----------
    %������������еĽ������ŵ����к����ź���ֵ����
    function stop = OutFun(x,optimValues,state)
        stop = false;   %ռλ����Ϊ�������ֵ����Ϊ��
           switch state     %���ݳ���ִ��״̬�б���Ҫִ�еĲ���
        %        case 'init'  %   fmincon������δ��ʼ
        %            hold on;
               case 'iter'  %   ĳ�ε�������
                    %++++++++++++++++++++++++++++
                    history.fval = [history.fval, optimValues.fval]; %��¼��ǰ�����Ľ������ź���ֵ
                    %++++++++++++++++++++++++++++
                    history.x = [history.x, x];  %��¼��ǰ�����Ľ������ŵ�
                    %nite = length(history.fval);   %��ǰ���������
                    %save('nite.dat','nite','-ascii'); 
%                     searchdir = [searchdir;...     %��¼ÿһ������������
%                        optimValues.searchdirection'];
%                 case 'done'  %fmincon��������
%                     hold off
               otherwise
           end
        end
    %Ŀ�꺯����Լ������
    %+++++++++++++++++++++++++++�����ݶ�+++++++++++++++++++++++++++++
    %Ŀ�꺯����Ŀ�꺯���ĵ���(������)
%     function [obj,gradobj] = pnonpobj (x)
%          obj = 2*sqrt(2) * x(1) + x(2);
%          gradobj = [2*sqrt(2);1];
%     end
%     %��ʽԼ���Ͳ���ʽԼ������(������)�Լ���ʽԼ���Ͳ���ʽԼ����������(����,��������Լ����������,���������Ż�����������)
%     function [conie,cone,gradconie,gradcone] = pnonpcon(x)
%         P = 2000; 
%         cone = [];
%         conie = [P * (sqrt(2) * x(1) + x(2)) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) - 2000,...
%                 -1500 + P * x(1) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2))];
%         gradcone = [];
%         gradconie = [P * (sqrt(2) * (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) - (2 * sqrt(2) * x(1) + 2 * x(2)) * (sqrt(2) * x(1) + x(2))) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) ^ 2,...
%                     P * ((sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) -  2 * x(1) * (sqrt(2) * x(1) + x(2))) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) ^ 2;...
%                     P * ((sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) - (2 * sqrt(2) * x(1) + 2 * x(2)) * x(1)) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) ^ 2,...
%                     - P * (2 * x(1) ^ 2) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) ^ 2];
%     end
    %+++++++++++++++++++++++++++����ݶ�+++++++++++++++++++++++++++++
    %Ŀ�꺯����Ŀ�꺯���ĵ���(������)
    % function [obj] = pnonpobj (x)
    %     obj = 2*sqrt(2) * x(1) + x(2);
    % end
    %��ʽԼ���Ͳ���ʽԼ������(������)�Լ���ʽԼ���Ͳ���ʽԼ����������(����,��������Լ����������,���������Ż�����������)
    % function [conie,cone] = pnonpcon(x)
    %    P = 2000; 
    %    cone = [];
    %    conie = [P * (sqrt(2) * x(1) + x(2)) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) - 2000,...
    %            -1500 + P * x(1) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2))];
    % end
    %++++++++++++++++++++++++++����Ԫ�����ṩԼ������+++++++++++++++++++++++++++
    %Ŀ�꺯����Ŀ�꺯���ĵ���(������)
    function [obj] = pnonpobj (x)
        obj = sqrt(2) * x(1) + x(2) + sqrt(2) * x(3);
    end
    %��ʽԼ���Ͳ���ʽԼ������(������)�Լ���ʽԼ���Ͳ���ʽԼ����������(����,��������Լ����������,���������Ż�����������)
    function [conie,cone] = pnonpcon(x)
        save('F:\WorkPath\ANSYS\SOP\cufile.dat','x','-ascii');
        cucommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Con_fun_update.py';
        system(cucommand);   %����python�ű����ύ����Ԫ����
        flg = 1;    %�������Ԫ�����Ƿ����н���
        while flg==1
            pause(1);
            if exist('F:\WorkPath\ANSYS\SOP\elemaxisstress.dat','file') == 2 %������Ԫ����ļ����ڣ�����Ϊ����Ԫ�������
                flg = 0;
            else
                if exist('F:\WorkPath\ANSYS\SOP\file.err','file') == 2
                    errcell = FileRead('F:\WorkPath\ANSYS\SOP\file.err',1000);  %��err�ļ��е��ַ������д���Ԫ������
                    for ii = 1 : 1 : length(errcell)
                        errloc = strfind(errcell{ii},'ERROR');  %���err�ļ����Ƿ���ڴ�����Ϣ
                        if ~isempty(errloc)     %��err�ļ��д��ڴ�����Ϣ������command windows����ʾ������Ϣ
                            disp('������ϢΪ��');
                            for kk = 1 : 1 : length(errcell)
                                disp(errcell{kk});
                            end
                            error('����Ԫ��������!');
                        end
                    end
                end
            end
        end
        %��ȡ����Ԫ������
        %axisf = load('F:\WorkPath\ANSYS\SOP\elemaxisforce.dat');
        %conie = [axisf(1)/ x(1) - 2000,axisf(2) / x(2) - 2000,-1500 - axisf(3) / x(3)];
        axst = load('F:\WorkPath\ANSYS\SOP\elemaxisstress.dat');
        conie = [axst(1) - 2000,axst(2) - 2000,-1500 - axst(3)];
        cone = []; 
        %����python�ű��������ļ���
        crcommand = 'Python F:\WorkPath\MATLAB\structural_optimization\RmTmpFile.py';   
        system(crcommand);
    end

end