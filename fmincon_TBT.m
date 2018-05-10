clear;clc;close all;
[history,timeconsumption,output] = runfmincon;
%clc;    %清屏，输出迭代信息
disp('=========================================');
disp('             程序执行信息：');
disp('-----------------------------------------')
disp(['总迭代次数：           ',num2str(output.iterations)]);
disp(['目标函数执行次数：     ',num2str(output.funcCount)]);
disp(['程序总运行时间：       ',num2str(timeconsumption),'  s']);
disp('最优点：');
disp('++++++++++++++++++++++');
disp(history.x(:,end));
disp('++++++++++++++++++++++');
disp(['最优函数值：            ',num2str(history.fval(end))]);
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
%将这一部分包含在函数内只是为了输出fmincon产生的近似最优点序列    
    x0 = [0.5;0.5;0.5];%初始点
    lb = [0.001;0.001;0.001];
    ub = [Inf;Inf;Inf];
    history.x = [];
    %+++++++++++++++++++++
    %当通过output function输出最优函数值序列或搜索方向序列时，需要将以下两条代码取消注释
    history.fval = [];
    %searchdir=[];
    %++++++++++++++++差分梯度++++++++++++++++
    options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed');  %fmincon函数的配置
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed');  %fmincon函数的配置
    %++++++++++++++++解析梯度++++++++++++++++
    %sqp算法
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed',...
    %    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);  %fmincon函数的配置
    %信赖域反射算法
    % options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter-detailed',...
    %     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);  %fmincon函数的配置
    %++++++++++++++++++++++++++++++++++++
    tstart = tic;   %开始计时
    %[x,fval,exitflag,output,lambda] = fmincon(@pnonpobj,x0,[],[],[],[],lb,ub,@pnonpcon,options);
    [~,~,exitflag,output,~] = fmincon(@pnonpobj,x0,[],[],[],[],lb,ub,@pnonpcon,options);
    if exitflag ~= 1
        disp(['exitflag =  ',num2str(exitflag)]);
        error('fmincon出错,请结合exitflag的值查找原因!');
    end
    tend = toc(tstart);     %计时结束
    timeconsumption = tend;
%========================================================================================
    %局部函数
%-----------
    %输出迭代过程中的近似最优点序列和最优函数值序列
    function stop = OutFun(x,optimValues,state)
        stop = false;   %占位，因为函数输出值不能为空
           switch state     %根据程序执行状态判别需要执行的操作
        %        case 'init'  %   fmincon迭代尚未开始
        %            hold on;
               case 'iter'  %   某次迭代结束
                    %++++++++++++++++++++++++++++
                    history.fval = [history.fval, optimValues.fval]; %记录当前迭代的近似最优函数值
                    %++++++++++++++++++++++++++++
                    history.x = [history.x, x];  %记录当前迭代的近似最优点
                    %nite = length(history.fval);   %当前迭代步编号
                    %save('nite.dat','nite','-ascii'); 
%                     searchdir = [searchdir;...     %记录每一步的搜索方向
%                        optimValues.searchdirection'];
%                 case 'done'  %fmincon迭代结束
%                     hold off
               otherwise
           end
        end
    %目标函数和约束函数
    %+++++++++++++++++++++++++++解析梯度+++++++++++++++++++++++++++++
    %目标函数和目标函数的导数(列向量)
%     function [obj,gradobj] = pnonpobj (x)
%          obj = 2*sqrt(2) * x(1) + x(2);
%          gradobj = [2*sqrt(2);1];
%     end
%     %等式约束和不等式约束函数(行向量)以及等式约束和不等式约束函数导数(矩阵,行数等于约束函数数量,列数等于优化变量的数量)
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
    %+++++++++++++++++++++++++++差分梯度+++++++++++++++++++++++++++++
    %目标函数和目标函数的导数(列向量)
    % function [obj] = pnonpobj (x)
    %     obj = 2*sqrt(2) * x(1) + x(2);
    % end
    %等式约束和不等式约束函数(行向量)以及等式约束和不等式约束函数导数(矩阵,行数等于约束函数数量,列数等于优化变量的数量)
    % function [conie,cone] = pnonpcon(x)
    %    P = 2000; 
    %    cone = [];
    %    conie = [P * (sqrt(2) * x(1) + x(2)) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2)) - 2000,...
    %            -1500 + P * x(1) / (sqrt(2) * x(1) ^ 2 + 2 * x(1) * x(2))];
    % end
    %++++++++++++++++++++++++++有限元计算提供约束函数+++++++++++++++++++++++++++
    %目标函数和目标函数的导数(列向量)
    function [obj] = pnonpobj (x)
        obj = sqrt(2) * x(1) + x(2) + sqrt(2) * x(3);
    end
    %等式约束和不等式约束函数(行向量)以及等式约束和不等式约束函数导数(矩阵,行数等于约束函数数量,列数等于优化变量的数量)
    function [conie,cone] = pnonpcon(x)
        save('F:\WorkPath\ANSYS\SOP\cufile.dat','x','-ascii');
        cucommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Con_fun_update.py';
        system(cucommand);   %调用python脚本，提交有限元任务
        flg = 1;    %检查有限元程序是否运行结束
        while flg==1
            pause(1);
            if exist('F:\WorkPath\ANSYS\SOP\elemaxisstress.dat','file') == 2 %若有限元结果文件存在，则认为有限元计算结束
                flg = 0;
            else
                if exist('F:\WorkPath\ANSYS\SOP\file.err','file') == 2
                    errcell = FileRead('F:\WorkPath\ANSYS\SOP\file.err',1000);  %将err文件中的字符串按行存入元胞数组
                    for ii = 1 : 1 : length(errcell)
                        errloc = strfind(errcell{ii},'ERROR');  %检查err文件中是否存在错误信息
                        if ~isempty(errloc)     %若err文件中存在错误信息，则在command windows中显示错误信息
                            disp('错误信息为：');
                            for kk = 1 : 1 : length(errcell)
                                disp(errcell{kk});
                            end
                            error('有限元分析出错!');
                        end
                    end
                end
            end
        end
        %读取有限元计算结果
        %axisf = load('F:\WorkPath\ANSYS\SOP\elemaxisforce.dat');
        %conie = [axisf(1)/ x(1) - 2000,axisf(2) / x(2) - 2000,-1500 - axisf(3) / x(3)];
        axst = load('F:\WorkPath\ANSYS\SOP\elemaxisstress.dat');
        conie = [axst(1) - 2000,axst(2) - 2000,-1500 - axst(3)];
        cone = []; 
        %调用python脚本，更新文件夹
        crcommand = 'Python F:\WorkPath\MATLAB\structural_optimization\RmTmpFile.py';   
        system(crcommand);
    end

end