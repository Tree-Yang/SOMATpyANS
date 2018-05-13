clear;clc;close all;
[history,timeconsumption,output] = runfmincon;
%clc;    %���������������Ϣ
if exist('sop_old.log','file') == 2
    delete('sop_old.log');   %ɾ��sop_old.log
end
if exist('sop.log','file') == 2
    movefile('sop.log','sop_old.log');   %��������־�ļ������Ƚ���������
end
diary('sop.log');    % ��־�ļ�Ŀ¼������
diary on;   %��ʼ�����־
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
diary off;  %���������־
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
legend('\itA_{\rm1}','\itA_{\rm2}','Location','SouthWest');
%plot(history.x(3,:),'-.d','LineWidth',1.5);
%legend('\itA_{\rm1}','\itA_{\rm2}','\itA_{\rm3}','Location','SouthWest');
grid on;
xlabel('Iteration');
ylabel('Section Area');
set(gca,'FontName','Euclid','FontSize',18);
%===========================================================================================
function [history,timeconsumption,output] = runfmincon
%����һ���ְ����ں�����ֻ��Ϊ�����fmincon�����Ľ������ŵ����� 
    x0 = [100;100];%��ʼ��
    lb = [0.001;0.001];
    ub = [500;500];
    history.x = [];
    %+++++++++++++++++++++
    %��ͨ��output function������ź���ֵ���л�������������ʱ����Ҫ��������������ȡ��ע��
    history.fval = [];
    %searchdir=[];
    %++++++++++++++++����ݶ�++++++++++++++++
    options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed','FiniteDifferenceType','central',...
        'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-4,'StepTolerance',1e-6);
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed','FiniteDifferenceType','central');%fmincon����������
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed','SpecifyObjectiveGradient',true);  %fmincon����������
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
    %+++++++++++++++++++++++++++����ݶ�+++++++++++++++++++++++++++++
    %Ŀ�꺯����Ŀ�꺯���ĵ���(������)
%     function [obj] = pnonpobj (x)
%          obj = (707 + 1000) * x(1) + (1414 + 1000) * x(2);
%          %gradobj = [1707;2414];
%     end
    %��ʽԼ���Ͳ���ʽԼ������(������)�Լ���ʽԼ���Ͳ���ʽԼ����������(����,��������Լ����������,���������Ż�����������)
%     function [conie,cone] = pnonpcon(x)
%         %xstr = 'x1';
%         %save('F:\WorkPath\MATLAB\cufile.dat',xstr,'-ascii');
%         %x = load('F:\WorkPath\MATLAB\cufile.dat');
%         l = [707, 1000, 1414, 1000];
%         NP = [-7500,-2500*sqrt(2),2500,5000*sqrt(2)];
%         N1 = [-2,-sqrt(2),1,sqrt(2)];
%         str = (NP ./ [x(1),x(1),x(2),x(2)]);
%         xx = [x(1),x(1),x(2),x(2)];
%         nd = 0;
%         for ii = 1:1:4
%             nd = nd + 1/(2.07e5) * N1(ii) * NP(ii) * l(ii) / xx(ii);
%         end
%         cone = [];
%         conie = [abs(str) - 100,abs(nd) - 1.5];
%     end
    %++++++++++++++++++++++++++����Ԫ�����ṩԼ������+++++++++++++++++++++++++++
    %Ŀ�꺯����Ŀ�꺯���ĵ���(������)
    function [obj] = pnonpobj (x)
         obj = (707 + 1000) * x(1) + (1414 + 1000) * x(2);
    end
    % %��ʽԼ���Ͳ���ʽԼ������(������)�Լ���ʽԼ���Ͳ���ʽԼ����������(����,��������Լ����������,���������Ż�����������)
    function [conie,cone] = pnonpcon(x)
        %clear classes;  %���MATLAB����python�����Ļ��棬�Ա�python�����ĸ���
        %����Ԫ����·��
        FEA_Path = 'F:\WorkPath\ANSYS\SOP\';
        %-----------------20180511 updated----------------------
        %����һ�γ���ʹ���ı���д�ķ�ʽ�����ݴ�matlab���ݵ�python�У����Ǵ��ںܴ�����⡣
        %��ʹ���ڽ���������������ļ���д�������ݵķ�ʽҲ�ᵼ�²�����
        %�����У���ֵ�����д�����ٶ�ȡ�õ�����ֵ֮��һ����1e-6��������
        %����Ϊ�����С�Ĳ��ᵼ���㷨�����һ�£�����û���ҵ�ԭ��
        %xstr = 'x';
        %save('F:\WorkPath\ANSYS\SOP\cufile.dat',xstr,'-ascii');
        %cucommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Con_fun_update.py';
        %system(cucommand);   %����python�ű����ύ����Ԫ����
        %--------------------By YJS------------------------------
        %--------------------20180512----------------------------
        %MATLAB�е���python�����ķ���
        %�ڽϸ߰汾��MATLAB�п���ֱ�ӵ���python��������Ҫ����׼�����ǽ�pythonĿ¼���뵽MATLAB��·����
        %py.**�ͱ�ʾ����python����
        py.ANSYS_mac_update.ANSYSmacupdate(x(1),x(1),x(2),x(2));
        %--------------------By YJS------------------------------
        system('Job_Submit.bat');
        flg = 1;    %�������Ԫ�����Ƿ����н���
        while flg==1
            pause(1);
            if exist([FEA_Path,'elemaxisstress.dat'],'file') == 2 %������Ԫ����ļ����ڣ�����Ϊ����Ԫ�������
                flg = 0;
            else
                if exist([FEA_Path,'Four_Bar_Truss.err'],'file') == 2
                    errcell = FileRead([FEA_Path,'Four_Bar_Truss.err'],1000);  %��err�ļ��е��ַ������д���Ԫ������
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
        axst = load([FEA_Path,'elemaxisstress.dat']);
        nd = load([FEA_Path,'nodedisp.dat']);
        conie = [abs(nd(3,2)) - 1.5,abs(axst)' - 100];
        cone = [];
        %����python�ű��������ļ���
        %----------------20180509------------------
        %�����������������п�����Ҫ�����������е�'Python'ȥ��
        %��Ϊ����һ���ű����Բ��ܰ���֮ǰ���ú����ķ�����ʵ��
        %����Ҳ���Գ��Խ���python�ű���дΪ������ʽ�ٵ���
        %+++++++++script++++++++++
        %crcommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Rm_tmp_file.py';   
        %system(crcommand);
        %+++++++++function++++++++
        py.Rm_tmp_file.RmTmpFile('F:\\WorkPath\\ANSYS\\');
        %------------------By YJS------------------
    end

end