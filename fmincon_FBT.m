clear;clc;close all;
[history,timeconsumption,output] = runfmincon;
%clc;    %清屏，输出迭代信息
if exist('sop_old.log','file') == 2
    delete('sop_old.log');   %删除sop_old.log
end
if exist('sop.log','file') == 2
    movefile('sop.log','sop_old.log');   %若存在日志文件，首先将其重命名
end
diary('sop.log');    % 日志文件目录及名称
diary on;   %开始输出日志
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
diary off;  %结束输出日志
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
%将这一部分包含在函数内只是为了输出fmincon产生的近似最优点序列 
    x0 = [100;100];%初始点
    lb = [0.001;0.001];
    ub = [500;500];
    history.x = [];
    %+++++++++++++++++++++
    %当通过output function输出最优函数值序列或搜索方向序列时，需要将以下两条代码取消注释
    history.fval = [];
    %searchdir=[];
    %++++++++++++++++差分梯度++++++++++++++++
    options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Algorithm','sqp','Display','iter-detailed','FiniteDifferenceType','central',...
        'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-4,'StepTolerance',1e-6);
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed','FiniteDifferenceType','central');%fmincon函数的配置
    %options = optimoptions(@fmincon,'OutputFcn',@OutFun,'Display','iter-detailed','SpecifyObjectiveGradient',true);  %fmincon函数的配置
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
    %+++++++++++++++++++++++++++差分梯度+++++++++++++++++++++++++++++
    %目标函数和目标函数的导数(列向量)
%     function [obj] = pnonpobj (x)
%          obj = (707 + 1000) * x(1) + (1414 + 1000) * x(2);
%          %gradobj = [1707;2414];
%     end
    %等式约束和不等式约束函数(行向量)以及等式约束和不等式约束函数导数(矩阵,行数等于约束函数数量,列数等于优化变量的数量)
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
    %++++++++++++++++++++++++++有限元计算提供约束函数+++++++++++++++++++++++++++
    %目标函数和目标函数的导数(列向量)
    function [obj] = pnonpobj (x)
         obj = (707 + 1000) * x(1) + (1414 + 1000) * x(2);
    end
    % %等式约束和不等式约束函数(行向量)以及等式约束和不等式约束函数导数(矩阵,行数等于约束函数数量,列数等于优化变量的数量)
    function [conie,cone] = pnonpcon(x)
        %clear classes;  %清除MATLAB调用python产生的缓存，以便python函数的更新
        %有限元分析路径
        FEA_Path = 'F:\WorkPath\ANSYS\SOP\';
        %-----------------20180511 updated----------------------
        %以下一段程序使用文本读写的方式将数据从matlab传递到python中，但是存在很大的问题。
        %即使对于解析函数的情况，文件读写传递数据的方式也会导致不收敛
        %测试中，数值本身和写出后再读取得到的数值之差一般在1e-6数量级。
        %至于为何如此小的差别会导致算法结果不一致，尚且没有找到原因。
        %xstr = 'x';
        %save('F:\WorkPath\ANSYS\SOP\cufile.dat',xstr,'-ascii');
        %cucommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Con_fun_update.py';
        %system(cucommand);   %调用python脚本，提交有限元任务
        %--------------------By YJS------------------------------
        %--------------------20180512----------------------------
        %MATLAB中调用python函数的方法
        %在较高版本的MATLAB中可以直接调用python函数，需要做的准备就是将python目录加入到MATLAB的路径中
        %py.**就表示调用python函数
        py.ANSYS_mac_update.ANSYSmacupdate(x(1),x(1),x(2),x(2));
        %--------------------By YJS------------------------------
        system('Job_Submit.bat');
        flg = 1;    %检查有限元程序是否运行结束
        while flg==1
            pause(1);
            if exist([FEA_Path,'elemaxisstress.dat'],'file') == 2 %若有限元结果文件存在，则认为有限元计算结束
                flg = 0;
            else
                if exist([FEA_Path,'Four_Bar_Truss.err'],'file') == 2
                    errcell = FileRead([FEA_Path,'Four_Bar_Truss.err'],1000);  %将err文件中的字符串按行存入元胞数组
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
        axst = load([FEA_Path,'elemaxisstress.dat']);
        nd = load([FEA_Path,'nodedisp.dat']);
        conie = [abs(nd(3,2)) - 1.5,abs(axst)' - 100];
        cone = [];
        %调用python脚本，更新文件夹
        %----------------20180509------------------
        %在其他电脑上运行有可能需要将下面命令中的'Python'去掉
        %因为这是一个脚本所以不能按照之前调用函数的方法来实现
        %或者也可以尝试将该python脚本改写为函数形式再调用
        %+++++++++script++++++++++
        %crcommand = 'Python F:\WorkPath\MATLAB\structural_optimization\Rm_tmp_file.py';   
        %system(crcommand);
        %+++++++++function++++++++
        py.Rm_tmp_file.RmTmpFile('F:\\WorkPath\\ANSYS\\');
        %------------------By YJS------------------
    end

end