file = open('time_small.txt','r')
file2 = open('friend_small.txt','r')
file3 = open('data_small.txt','w')
data1 = file.readlines()
data2 = file2.readlines()
for i in range(len(data1)):
    file3.write(data1[i][:-1]+' '+data2[i][:-1]+'\n')
file.close()
file2.close()
file3.close()
