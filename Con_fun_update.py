#更新约束函数
import os
#从文件读取输入数据
fcu_data = []
with open('F:\WorkPath\ANSYS\SOP\cufile.dat', mode = 'r', encoding = 'utf8') as fcu:
    fcu_lines = fcu.readlines()  
for line in fcu_lines:  
    lineArr = line.strip()
    fcu_data.append(lineArr)
secarea1 = float(fcu_data[0])
secarea2 = float(fcu_data[1])
secarea3 = float(fcu_data[2])

#----------20180510-----------
#不需要迭代次数，直接在有限元工作目录下创建新文件
#by  YJS
# fni_data = []
# with open('nite.dat', mode = 'r', encoding = 'utf8') as fni:
#     fni_lines = fni.readlines()  
# for line in fni_lines:  
#     lineArr1 = line.strip()
#     fni_data.append(lineArr1)

# nite = int(float(fni_data[0]))

#更新有限元模型
from ANSYS_mac_update import ANSYSmacupdate
#from Job_Submit_update import JobSubmitUpdate
if __name__ == "__main__":
    ANSYSmacupdate(secarea1,secarea2,secarea3)
    #JobSubmitUpdate(nite)
#运行有限元模型
os.system('F:\\WorkPath\\MATLAB\\structural_optimization\\Job_Submit.bat')
