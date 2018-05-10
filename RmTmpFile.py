#import shutil
import os
#if os.path.exists('F:\\WorkPath\\ANSYS\\SOP\\elemaxisforce.dat'):
#    shutil.move('F:\\WorkPath\\ANSYS\\SOP\\elemaxisforce.dat','F:\\WorkPath\\MATLAB\\structural_optimization\
#\\elemaxisforce.dat') 
#shutil.rmtree('F:\\xxx')  
path = 'F:\\WorkPath\\ANSYS\\'
count = 0
for fn in os.listdir(path): #fn 表示的是文件名,计算当前目录下的文件夹数量
        count = count + 1
#print(count)
ncurdir = count - 1
diroldname = 'SOP' + str(ncurdir).zfill(4)
diroldpath = path + diroldname
os.rename('F:\\WorkPath\\ANSYS\\SOP',diroldpath)        #重命名使用过的SOP文件夹
os.mkdir('F:\\WorkPath\\ANSYS\\SOP')    #创建新的SOP文件夹