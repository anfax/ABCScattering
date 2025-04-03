import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # 导入3D绘图工具
from collections import defaultdict

# 设置目录
dir = './output'
fil_e = os.path.join(dir, 'energy.dat')
dat_e = np.loadtxt(fil_e)
energy = dat_e[:, 1]

 
files = os.listdir(dir)
files = [f for f in os.listdir(dir) if f.startswith('S-matrix') and f.endswith('.dat')]
# files.sort(key=lambda x: int(x.split('-')[1].split('.')[0]))
label_state = []
for ii in range(1,6): 
    ee = [] 
    ss = [] 
    # 遍历文件并绘制图形
    # for file in files:
    for i in range(1,100):

        # 拼接文件路
        file = f'S-matrix{i}.dat'
        fil = os.path.join(dir, file)
        
        # 加载数据，跳过第一行（假设是标题行）
        dat = np.loadtxt(fil)
       
            
        # 提取数据列
        a = dat[:, 0]
        v = dat[:, 1]
        j = dat[:, 2]
        k = dat[:, 3]
        aa = dat[:,4]
        
        vv = dat[:, 5]
        jj = dat[:, 6]
        kk = dat[:, 7]
        r = dat[:, 8]
        i_val = dat[:, 9]  # 避免与变量i冲突，改名为i_val
        sq = dat[:, 10]
        ee.append(e[i])
        ss.append(sq[ii])
        
    # 显示图形
    
    plt.plot(ee,ss,label=f'$|{int(a[ii])},{int(v[ii])},{int(j[ii])},{int(k[ii])}\\rangle -> |{int(aa[ii])},{int(vv[ii])},{int(jj[ii])},{int(kk[ii])}\\rangle$')
    # plt.tight_layout()
    
plt.legend()
plt.savefig(f'fig/smatrixs.jpg',dpi=666)
plt.close()
        
