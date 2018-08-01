import pandas as pd
#data = pd.read_csv("index.csv")
#file = open("index.txt",'w')
#for i in range(len(data)):
#	file.write(str(data.iloc[i]['post_time_day'])+' ')
#	file.write(str(int(data.iloc[i]['start_ind']))+' ')
#	file.write(str(int(data.iloc[i]['end_ind']))+'\n')
#file.close()
data = pd.read_csv("data.csv")
file = open("data.txt",'w')
for i in range(len(data)):
	file.write(str(int(data.iloc[i]['relative_time_second']))+' ')
	file.write(str(int(data.iloc[i]['number_of_followers']))+'\n')
file.close()

