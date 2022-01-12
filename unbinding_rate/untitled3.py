# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:20:51 2022

@author: stefa
"""

import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


folder =r"C:\Users\stefa\OneDrive - University of Wollongong\phi29Paper\data\lambda\211221\norm_integration/"
files = os.listdir(folder)

data = pd.read_csv(folder+files[0])


data['seconds'] = data['slice'] * 0.2


data = data[data['seconds']>0]

grouped = data.groupby('trajectory')

traj_numbers = data['trajectory'].unique()


n1=543
n2=27
n3=316

traj1 = grouped.get_group(n1)
traj2 = grouped.get_group(n2)
traj3 = grouped.get_group(n3)


#print(traj1.keys())

fig1,ax1 = plt.subplots()
ax1.plot(traj1['seconds'],traj1['nts'],label = 'traj {}'.format(n1),color = 'black')

ax1.legend()
ax1.set_xlim([25,450])
#fig1.savefig('Lambda_traj{}.svg'.format(n1))


# fig2,ax2 = plt.subplots()
# ax2.plot(traj2['seconds'],traj2['nts'],label = 'traj {}'.format(n2),color = 'black')
# ax2.legend()
# ax2.set_xlim([50,400])
# #fig2.savefig('Lambda_traj{}.svg'.format(n2))


# fig3,ax3 = plt.subplots()
# ax3.plot(traj3['seconds'],traj3['nts'],label = 'traj {}'.format(n3),color = 'black')
# ax3.legend()
# ax3.set_xlim([50,400])
#fig3.savefig('Lambda_traj{}.svg'.format(n3))



plt.show()

choices = [269,79,105,22,220,316,152,2,120]