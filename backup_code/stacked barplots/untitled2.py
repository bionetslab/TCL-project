# libraries
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
 
# set width of bars
barWidth = 0.21
 
# set heights of bars
NB_Means = [84.76, 81.68, 84.45, 86.54, 87.01]
SVM_Means = [83.35, 82.21, 82.90, 83.17, 87.42]
LDA_Means = [84.86, 84.39, 84.22, 83.95, 89.39]
DT_Means = [78.77, 80.01, 77.49, 84.16, 84.66]

NB_Std = np.array([0.24*2, 0.39*2, 0.25*2, 0.30*2, 0.30*2])
SVM_Std = np.array([0.43*2, 0.43*2, 0.31*2, 0.40*2, 0.33*2])
LDA_Std = np.array([0.26*2, 0.27*2, 0.09*2, 0.22*2, 0.26*2])
DT_Std = np.array([0.89*2, 1.16*2, 1.10*2, 0.89*2, 0.86*2])

# NB_Std = np.array([4, 3, 4, 1, 5])
# SVM_Std = np.array([3, 5, 2, 3, 3])
# LDA_Std = np.array([4, 3, 4, 1, 5])
# DT_Std = np.array([3, 5, 2, 3, 3])
 
# Set position of bar on X axis
r1 = np.arange(len(NB_Means))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]
 
# Make the plot
rects1=plt.bar(r1, NB_Means, width=barWidth, yerr=NB_Std, label='NB')
rects2=plt.bar(r2, SVM_Means, width=barWidth, yerr=SVM_Std, label='SVM')
rects3=plt.bar(r3, LDA_Means, width=barWidth, yerr=LDA_Std, label='LDA')
rects4=plt.bar(r4, DT_Means, width=barWidth, yerr=DT_Std, label='DT')
 
# # Add xticks on the middle of the group bars
# # sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
# # Customize the major grid
# plt.grid(which='major', linestyle='-', linewidth='0.5', color='red')
# # Customize the minor grid
# plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

plt.grid(axis='y', which='major', color='k', linestyle='-')
plt.grid(axis='y', which='minor', color='r', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.rcParams['axes.axisbelow'] = True

plt.axis([None, None, 0, 100])
plt.title('Normal vs Cancer')
plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(NB_Means))], ['AlexNet', 'GoogleNet', 'InceptionV3', 'ResNet50', 'XceptionNet'])


# plt.show()

# plt.bar_label(rects1, padding=3)
# plt.bar_label(rects2, padding=3)
# plt.bar_label(rects3, padding=3)
# plt.bar_label(rects4, padding=3)
 
# Create legend & Show graphic
plt.legend()
plt.show()






















