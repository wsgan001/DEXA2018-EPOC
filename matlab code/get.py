file = open('data.txt','r')
t = 0
for i in file.readlines():
    t += 1
print t
file.close()
