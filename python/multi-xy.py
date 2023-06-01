import matplotlib.pyplot as plt

# Read data from file1.dat to file4.dat
with open('fisbar.dat', 'r') as file1:
    next(file1)  # Skip first line (comment)
    data1 = [[float(x) for x in line.strip().split()] for line in file1]

with open('fusbar.dat', 'r') as file2:
    next(file2)  # Skip first line (comment)
    data2 = [[float(x) for x in line.strip().split()] for line in file2]

with open('xy3.dat', 'r') as file3:
    next(file3)  # Skip first line (comment)
    data3 = [[float(x) for x in line.strip().split()] for line in file3]

with open('xy4.dat', 'r') as file4:
    next(file4)  # Skip first line (comment)
    data4 = [[float(x) for x in line.strip().split()] for line in file4]

with open('xy5.dat', 'r') as file5:
    next(file5)  # Skip first line (comment)
    data5 = [[float(x) for x in line.strip().split()] for line in file5]

with open('xy6.dat', 'r') as file6:
    next(file6)  # Skip first line (comment)
    data6 = [[float(x) for x in line.strip().split()] for line in file6]

# Create a 2x3 subplot grid
fig, axs = plt.subplots(2, 3)

# Plot data on each subplot
axs[0, 0].plot([row[0] for row in data1], [row[1] for row in data1])
axs[0, 0].set_title('fisbar')
#axs[0, 0].set_xlim([5, 30])

axs[0, 1].plot([row[0] for row in data2], [row[1] for row in data2])
axs[0, 1].set_title('fusbar')
#axs[0, 1].set_xlim([5, 30])

axs[0, 2].plot([row[0] for row in data3], [row[1:] for row in data3])
axs[0, 2].set_title('Data 3')

axs[1, 0].plot([row[0] for row in data4], [row[1:] for row in data4])
axs[1, 0].set_title('Data 4')

axs[1, 1].plot([row[0] for row in data5], [row[1:] for row in data5])
axs[1, 1].set_title('Data 5')

axs[1, 2].plot([row[0] for row in data6], [row[1:] for row in data6])
axs[1, 2].set_title('Data 6')

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.5)

# Show the plot
plt.savefig('multifig.pdf', format='pdf')
plt.show()
