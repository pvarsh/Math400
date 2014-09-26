from math import * 

n = 50
mylist = []
for i in range(n):
    mylist.append(i)

print mylist


mylist1 = [sin(i) + cos(2*i) for i in range(n)]
print mylist1
