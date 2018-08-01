file = open('index.txt','r')
file2 = open('data.txt','r')
file3 = open('index_small.txt','w')
file4 = open('time_small.txt','w')
file5 = open('friend_small.txt','w')
data = file2.readlines()
for i in file.readlines()[0:10000]:
    file3.write(i)
    temp = i.split(' ')
    for j in data[int(temp[1])-1:int(temp[2])]:
        temp2 = j.split(' ')
        file4.write(str(temp2[0])+'\n')
        file5.write(str(temp2[1][:-2])+'\n')
file.close()
file2.close()
file3.close()
file4.close()
file5.close()
    
    
