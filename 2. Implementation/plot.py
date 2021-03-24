import matplotlib.pyplot as plt


X, Y = [], []
for line in open('PFPASAP_RES.txt', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])

plt.plot(X, Y, '-o')
plt.title('plot for utilization ')
plt.xlabel('utilization')
plt.ylabel('percent')
plt.show()
