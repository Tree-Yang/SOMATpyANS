#该python函数用来更新ANSYS的mac文件
def ANSYSmacupdate(secarea1,secarea2,secarea3,secarea4):
    #import os
    #++++++++++++++++++++++++++++++++++++++++++
    #Updated by YJS 20180508
    #其实不需要给mac文件编号
    #nite_old = nite - 1
    #file_name_old = 'Tri_Bar_Truss' + str(nite_old).zfill(4) + '.mac'   #上一次迭代使用的mac文件
    #file_path_old = 'F:\\WorkPath\\MATLAB\\structural_optimization\\' + file_name_old
    #file_name = 'Tri_Bar_Truss' + str(nite).zfill(4) + '.mac'   #需要新生成的mac文件
    #file_path = 'F:\\WorkPath\\MATLAB\\structural_optimization\\' + file_name
    #if os.path.exists('Tri_Bar_Truss_old.mac'): #若该文件存在则删除该文件，否则将导致下一步的文件重命名失败
    #    os.remove('Tri_Bar_Truss_old.mac')
    #os.rename('Tri_Bar_Truss.mac','Tri_Bar_Truss_old.mac')  #重命名
    #file_name_old = 'Tri_Bar_Truss_old.mac'   #上一次迭代使用的mac文件
    file_name = 'FBT.mac'   #需要新生成的mac文件
    file_path_old = 'F:/WorkPath/MATLAB/structural_optimization/FBT/' + file_name
    file_path = 'F:/WorkPath/ANSYS/SOP/' + file_name
    #++++++++++++++++++++++++++++++++++++++\++++++
    #需要指定编码格式为utf8，否则可能会报错"UnicodeDecodeError".
    with open(file_path_old, mode = 'r', encoding = 'utf8') as f_r: #读入待修改内容
        f_r_lines = f_r.readlines()
    f_r_lines_list = []
    for line in f_r_lines:  #文本内容按行放入列表
        f_r_lines_list.append(line)
    #新内容
    #print('secarea1 = ' + str(secarea1) + '$secarea2 = ' + str(secarea2) + '$secarea3 = ' + str(secarea3))
    f_r_lines_list[14] = 'secarea1 = ' + str(secarea1) + '$secarea2 = ' + str(secarea2) + '$secarea3 = ' + str(secarea3)+'$secarea4 = ' + str(secarea4)+'\n'
    #修改文本内容，写入新文件
    #--------------------method 1--------------------
    #with open(file_path, mode = 'w', encoding = 'utf8') as f_w:
    #    for ii in range(len(f_r_lines_list)):
    #        f_w.write(f_r_lines_list[ii])
    #--------------------method 2--------------------
    with open(file_path, mode = 'w', encoding = 'utf8') as f_w:
        f_w.writelines(f_r_lines_list)
    return

#测试用
ANSYSmacupdate(1.059112482670865e+02,82.942677736246880,1.059112482670865e+02,82.942677736246880)