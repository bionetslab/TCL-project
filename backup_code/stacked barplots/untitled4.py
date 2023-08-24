# Load Matplotlib and data wrangling libraries.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load jobs dataset from Vega's dataset library.
from vega_datasets import data

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
 


# Use Seaborn's context settings to make fonts larger.
import seaborn as sns
sns.set_context('talk')

# Create a grouped bar chart, with job as the x-axis
# and gender as the variable we're grouping on so there
# are two bars per job.
fig, ax = plt.subplots(figsize=(12, 8))

# Our x-axis. We basically just want a list
# of numbers from zero with a value for each
# of our jobs.
# x = np.arange(len(df.job.unique()))

# Set position of bar on X axis
r1 = np.arange(len(NB_Means))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]

# Define bar width. We need this to offset the second bar.
# bar_width = 0.4

# Make the plot
rects1=ax.bar(r1, NB_Means, width=barWidth, yerr=NB_Std, label='NB', capsize=4)
rects2=ax.bar(r2, SVM_Means, width=barWidth, yerr=SVM_Std, label='SVM', capsize=4)
rects3=ax.bar(r3, LDA_Means, width=barWidth, yerr=LDA_Std, label='LDA', capsize=4)
rects4=ax.bar(r4, DT_Means, width=barWidth, yerr=DT_Std, label='DT', capsize=4)

# b1 = ax.bar(x, df.loc[df['sex'] == 'men', 'count'],
#             width=bar_width, label='Men')
# # Same thing, but offset the x.
# b2 = ax.bar(x + bar_width, df.loc[df['sex'] == 'women', 'count'],
#             width=bar_width, label='Women')

# # Fix the x-axes.
# ax.set_xticks(x + bar_width / 2)
# ax.set_xticklabels(df.job.unique())

# Add legend.
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.105),
          ncol=4, fancybox=True, shadow=True, fontsize=16)

plt.grid(axis='y', which='major', color='#BEBEBE', linestyle='-')
plt.grid(axis='y', which='minor', color='#BEBEBE', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.rcParams['axes.axisbelow'] = True

plt.axis([None, None, 0, 100])
plt.title('Normal vs Cancer\n\n', fontweight='bold')
plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(NB_Means))], ['AlexNet', 'GoogleNet', 'InceptionV3', 'ResNet50', 'XceptionNet'])



for tick in ax.xaxis.get_major_ticks():
    tick.tick1line.set_visible(True)
    tick.tick2line.set_visible(False)
    tick.label1.set_visible(True)
    tick.label2.set_visible(False)


for tick in ax.xaxis.get_minor_ticks():
    tick.tick1line.set_visible(False)





for bar in ax.patches:
  # The text annotation for each bar should be its height.
  bar_value = bar.get_height()
  # Format the text with commas to separate thousands. You can do
  # any type of formatting here though.
  text = f'   {bar_value:,}'
  # This will give the middle of each bar on the x-axis.
  text_x = bar.get_x() + bar.get_width() / 2
  # get_y() is where the bar starts so we add the height to it.
  text_y = bar.get_y() + bar_value
  # If we want the text to be the same color as the bar, we can
  # get the color like so:
  bar_color = bar.get_facecolor()
  # If you want a consistent color, you can just set it as a constant, e.g. #222222
  ax.text(text_x, text_y, text, ha='center', va='bottom', color=bar_color,
          size=12, rotation=90)

# # Axis styling.
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_color('#DDDDDD')
# ax.tick_params(bottom=False, left=False)
# ax.set_axisbelow(True)
# ax.yaxis.grid(True, color='#EEEEEE')
# ax.xaxis.grid(False)

# # Add axis and chart labels.
# ax.set_xlabel('Job', labelpad=15)
# ax.set_ylabel('# Employed', labelpad=15)
# ax.set_title('Employed Workers by Gender for Select Jobs', pad=15)

fig.tight_layout()