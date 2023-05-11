import matplotlib.pyplot as plt

# 读取数据文件
with open('./fusbar.dat', 'r') as f:
    data = f.readlines()
#with open('./fisbar_def.dat', 'r') as f2:
#    data2 = f2.readlines()

# 提取x和y坐标
x1, y1, y2 = [], [], []
for line in data[1:]:
    line = line.strip().split()
    x1.append(float(line[0]))
    y1.append(float(line[1]))
    y2.append(float(line[2]))

#x12, y12, y22 = [], [], []
#for line in data2[222:]:
#    line = line.strip().split()
#    x12.append(float(line[0]))
#    y12.append(float(line[1]))
#    y22.append(float(line[2]))

# 绘制折线图
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('R (fm)')
ax1.set_ylabel('E (MeV)')
ax1.plot(x1, y1, color=color, label='energy')
#ax.tick_params(axis='y', labelcolor=color)

#color2 = 'tab:blue'
#ax1.set_xlabel('R (fm)')
#ax1.set_ylabel('E (MeV)')
#ax1.plot(x12, y12, color=color2, label='energy2')
#ax1.tick_params(axis='y', labelcolor=color)

#使用ax1.twinx()来创建一个与ax1共享x轴但具有不同y轴的轴ax2
#ax2 = ax1.twinx()
#color = 'tab:blue'
#ax2.set_ylabel('DER', color=color)
#ax2.plot(x1, y2, color=color)
#ax2.tick_params(axis='y', labelcolor=color, label='der')

# 添加图例
lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2, loc='upper right')
ax1.legend(lines, labels, loc='upper right')


plt.title('α decay')

fig.tight_layout()
plt.savefig('fisbar.pdf', format='pdf')
plt.show()
