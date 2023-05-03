import matplotlib.pyplot as plt

# 读取数据文件
with open('./fort.7', 'r') as f:
    data = f.readlines()


# 提取x和y坐标
x1, y1, y2, y3 = [], [], [], []
for line in data[:]:
    line = line.strip().split()
    x1.append(float(line[0]))
    y1.append(float(line[1]))
    y2.append(float(line[2]))
    y3.append(float(line[3]))


# 绘制折线图
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('R (fm)')
ax1.set_ylabel('E (MeV)')
ax1.plot(x1, y1, color=color, label='energy')
ax1.plot(x1, y2, color='tab:blue', label='energy2')
ax1.plot(x1, y3, color='tab:green', label='energy3')
#ax1.tick_params(axis='y', labelcolor=color)



#使用ax1.twinx()来创建一个与ax1共享x轴但具有不同y轴的轴ax2
#ax2 = ax1.twinx()
#color = 'tab:blue'
#ax2.set_ylabel('DER', color=color)
#ax2.plot(x1, y2, color=color, label='energy2')
#ax2.plot(x1, y3, color='tab:green', label='energy3')
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
