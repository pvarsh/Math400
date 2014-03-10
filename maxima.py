x = [0, -1, 2, 3, -10, 7, 40, -50, 90, 10]
#print abs(max(x[0:-1], key = lambda x: abs(x)))

#print abs(max(x[0:-1], key = lambda x: abs(x)))
print x[0:4]

maxIndex = max(range(len(x)), key = lambda i: abs(x[i]))
print maxIndex
