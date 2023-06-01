import matplotlib.pyplot as plt
import numpy as np

# 创建一个2x3的点线图子图
fig, axs = plt.subplots(2, 3)

# 定义每个子图的x和y坐标范围
x1 = np.linspace(0, 10, 100)
y1 = np.sin(x1)
axs[0, 0].plot(x1, y1)
axs[0, 0].set_xlim([0, 10])
axs[0, 0].set_ylim([-1, 1])

x2 = np.linspace(0, 5, 50)
y2 = np.exp(x2)
axs[0, 1].plot(x2, y2)
axs[0, 1].set_xlim([0, 5])
axs[0, 1].set_ylim([0, 150])

x3 = np.linspace(-5, 5, 50)
y3 = x3 ** 2
axs[0, 2].plot(x3, y3)
axs[0, 2].set_xlim([-5, 5])
axs[0, 2].set_ylim([0, 25])

x4 = np.linspace(-10, 10, 200)
y4 = np.tan(x4)
axs[1, 0].plot(x4, y4)
axs[1, 0].set_xlim([-10, 10])
axs[1, 0].set_ylim([-10, 10])

x5 = np.linspace(0, 20, 100)
y5 = np.log(x5)
axs[1, 1].plot(x5, y5)
axs[1, 1].set_xlim([0, 20])
axs[1, 1].set_ylim([-5, 5])

x6 = np.linspace(-20, 20, 200)
y6 = np.sqrt(abs(x6))
axs[1, 2].plot(x6, y6)
axs[1, 2].set_xlim([-20, 20])
axs[1, 2].set_ylim([0, 5])

# 显示图形
plt.show()
