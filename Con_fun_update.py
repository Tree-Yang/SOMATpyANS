#更新约束函数
import os

#从文件读取输入数据
fcu_data = []
with open('cufile.dat', mode = 'r', encoding = 'utf8') as fcu:
    fcu_lines = fcu.readlines()  
for line in fcu_lines:  
    lineArr = line.strip()
    fcu_data.append(lineArr)
secarea1 = float(fcu_data[0])
secarea2 = float(fcu_data[1])
secarea3 = secarea1

#更新有限元模型
from ANSYS_mac_update import ANSYSmacupdate
if __name__ == "__main__":
    ANSYSmacupdate(secarea1,secarea2,secarea3)

#运行有限元模型
os.system('F:\\WorkPath\\MATLAB\\structural_optimization\\Job_Submit.bat')