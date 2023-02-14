import pandas as pd
import numpy as np
import datetime

n_dir = {"ULK1": 7, "KRAS": 238, "BCL2": 29, "MHL1": 8, 
         "PPIA": 152, "BCL2L1": 97, "TRAP1": 20, "ARID3A": 3, 
         "MAP2LC3B": 25, "FAK1": 35, "FKBP1A": 59, "MZF1": 2}

n_structure = sum(n_dir.values())

df = pd.DataFrame(columns=['cores','protein',
                       'start', 'total', "per_pdb"])

files = ["ARID3A_log.txt", "BCL2L1_log.txt", "BCL2_log.txt", "FAK1_log.txt",
         "FKBP1A_log.txt", "KRAS_log.txt", "MAP2LC3B_log.txt", "MHL1_log.txt",
         "MZF1_log.txt","PPIA_log.txt", "TRAP1_log.txt", "ULK1_log.txt"]

for i in range(1,5):
    with open(f'{i}core/log.txt') as f:
        lines = f.readlines()
    
    starttime = lines[6].split(": ")[1].split("\n")[0] 
    (h, m, s) = starttime.split(':')
    st = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))

    totaltime = lines[5].split(": ")[1].split("\n")[0]
    (h, m, s) = totaltime.split(':')
    t = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
    
    per_structure = (t - st)/n_structure
    
    list_of_values = [i, "all", starttime, totaltime, per_structure]
    
    df.loc[i] = np.array(list_of_values, dtype="object")
    
for i in files:
    j = i.split("_")[0]
    with open(f"commandline/{i}") as f:
        lines = f.readlines()
    
    starttime = lines[6].split(": ")[1].split("\n")[0] 
    (h, m, s) = starttime.split(':')
    st = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))

    totaltime = lines[5].split(": ")[1].split("\n")[0]
    (h, m, s) = totaltime.split(':')
    t = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
    
    per_structure = (t - st)/n_dir[j]

    list_of_values = [1, j, starttime, totaltime, per_structure]
    
    df.loc[i] = np.array(list_of_values, dtype="object")

df1 = df[5:]
timeList = list(df1.total)
mysum = datetime.timedelta()
for i in timeList:
    (h, m, s) = i.split(':')
    d = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
    mysum = mysum + d

total_time = mysum/60

per_structure = mysum/n_structure

list_of_values = [1, "all", "NA", total_time, per_structure]

df.loc["total"] = np.array(list_of_values, dtype="object")

df.to_csv("time.csv")
